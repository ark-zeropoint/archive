#!/usr/bin/env python3
"""
Taylor-Green Vortex Re=1600 benchmark on T^3 = [0, 2*pi]^3.

Three cases, same pseudospectral RK4 integrator:
  1. Unresolved baseline DNS          (nu = 1/1600, no model)
  2. Standard Smagorinsky LES         (Cs = 0.17)
  3. NS_TZT selective dissipation     (sealed Phase-3 operators)

Outputs
  tgv_metrics.csv
  tgv_dissipation_rate.png
  tgv_energy_spectrum.png
  tgv_max_vorticity.png

Run
  python tzt_tgv_benchmark.py
  python tzt_tgv_benchmark.py --n 32 --t-end 2.0     # smoke
  python tzt_tgv_benchmark.py --n 128 --dt 0.002     # fine (slow in pure NumPy)
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from numpy.fft import fftn, ifftn

# ---------------------------------------------------------------------------
# Grid and spectral operators
# ---------------------------------------------------------------------------

def make_grid(n: int):
    """Periodic [0, 2*pi]^3, full complex FFT layout."""
    x = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    k = np.fft.fftfreq(n, d=1.0 / n).astype(np.float64)  # integer wavenumbers
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K2 = KX * KX + KY * KY + KZ * KZ
    K = np.sqrt(K2)
    # 2/3 dealiasing mask
    kcut = n / 3.0
    mask = ((np.abs(KX) < kcut) & (np.abs(KY) < kcut) & (np.abs(KZ) < kcut)).astype(np.float64)
    return X, Y, Z, KX, KY, KZ, K2, K, mask


def project(uh, vh, wh, KX, KY, KZ, K2):
    """Leray projector. Zero mode is left untouched."""
    # k · u
    div = KX * uh + KY * vh + KZ * wh
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(K2 > 0.0, 1.0 / K2, 0.0)
    uh = uh - KX * div * inv
    vh = vh - KY * div * inv
    wh = wh - KZ * div * inv
    uh = uh.copy(); vh = vh.copy(); wh = wh.copy()
    uh[0, 0, 0] = 0.0
    vh[0, 0, 0] = 0.0
    wh[0, 0, 0] = 0.0
    return uh, vh, wh


def curl_hat(uh, vh, wh, KX, KY, KZ):
    """Fourier curl: i k × û."""
    i = 1j
    cx = i * (KY * wh - KZ * vh)
    cy = i * (KZ * uh - KX * wh)
    cz = i * (KX * vh - KY * uh)
    return cx, cy, cz


def grad_hat(fh, KX, KY, KZ):
    i = 1j
    return i * KX * fh, i * KY * fh, i * KZ * fh


def dealias(uh, vh, wh, mask):
    return uh * mask, vh * mask, wh * mask


def to_phys(fh):
    return ifftn(fh).real


def to_hat(f):
    return fftn(f)


def energy(uh, vh, wh, n):
    """E = 0.5 < |u|^2 >. fftn is unnormalized, so divide by n^3 twice via Parseval once."""
    scale = 1.0 / n**3
    e = 0.5 * scale * (np.sum(np.abs(uh) ** 2) + np.sum(np.abs(vh) ** 2) + np.sum(np.abs(wh) ** 2)) / n**3
    # Wait: Parseval for numpy fftn: sum(|u|^2) = sum(|û|^2) / n^3
    # so mean(|u|^2) = sum(|û|^2) / n^6
    # Correct:
    return e  # already mean because sum/n^3 / n^3 = sum/n^6, times 0.5. Yes.


def max_vorticity(uh, vh, wh, KX, KY, KZ):
    cx, cy, cz = curl_hat(uh, vh, wh, KX, KY, KZ)
    wx, wy, wz = to_phys(cx), to_phys(cy), to_phys(cz)
    return float(np.max(np.sqrt(wx * wx + wy * wy + wz * wz)))


def shell_spectrum(uh, vh, wh, K, n):
    """3D shell-averaged energy spectrum E(k), integer shells."""
    scale = 0.5 / n**6  # so sum_k E(k) ~= kinetic energy
    kk = np.rint(K).astype(np.int32)
    kmax = n // 2
    ek = np.zeros(kmax + 1, dtype=np.float64)
    spec = scale * (np.abs(uh) ** 2 + np.abs(vh) ** 2 + np.abs(wh) ** 2)
    # bin
    for s in range(0, kmax + 1):
        ek[s] = np.sum(spec[kk == s])
    return ek


# ---------------------------------------------------------------------------
# Initial condition: Taylor-Green
# ---------------------------------------------------------------------------

def taylor_green(X, Y, Z):
    u = np.sin(X) * np.cos(Y) * np.cos(Z)
    v = -np.cos(X) * np.sin(Y) * np.cos(Z)
    w = np.zeros_like(u)
    return u, v, w


# ---------------------------------------------------------------------------
# Model forces (physical -> Fourier, then projected by caller)
# ---------------------------------------------------------------------------

def phi_switch(xi, beta):
    """Smooth gate. Hard zero for xi < 0.8 (laminar lock)."""
    phi = 0.5 * (1.0 + np.tanh(beta * (xi - 1.2)))
    phi = np.where(xi < 0.8, 0.0, phi)
    return phi


def force_smagorinsky(uh, vh, wh, KX, KY, KZ, n, cs, mask):
    """Standard dissipative Smagorinsky.
    nu_t = (Cs * Delta)^2 * |S|,  F = div(2 nu_t S).
    User's tau_ij = 2 (Cs Delta)^2 |S| S_ij is the stress magnitude;
    the momentum force is +div(tau) with this sign so that energy decays.
    """
    delta = 2.0 * np.pi / n
    # gradients
    u, v, w = to_phys(uh), to_phys(vh), to_phys(wh)
    # spectral derivatives of physical fields via already-dealiased hats
    ux, uy, uz = [to_phys(g) for g in grad_hat(uh, KX, KY, KZ)]
    vx, vy, vz = [to_phys(g) for g in grad_hat(vh, KX, KY, KZ)]
    wx, wy, wz = [to_phys(g) for g in grad_hat(wh, KX, KY, KZ)]
    # strain
    sxx = ux
    syy = vy
    szz = wz
    sxy = 0.5 * (uy + vx)
    sxz = 0.5 * (uz + wx)
    syz = 0.5 * (vz + wy)
    s2 = sxx * sxx + syy * syy + szz * szz + 2.0 * (sxy * sxy + sxz * sxz + syz * syz)
    s_mag = np.sqrt(2.0 * s2 + 1e-30)
    nu_t = (cs * delta) ** 2 * s_mag
    # stress
    txx, tyy, tzz = 2.0 * nu_t * sxx, 2.0 * nu_t * syy, 2.0 * nu_t * szz
    txy, txz, tyz = 2.0 * nu_t * sxy, 2.0 * nu_t * sxz, 2.0 * nu_t * syz
    # divergence in Fourier
    def ddx(f):
        return to_hat(f) * (1j * KX)
    def ddy(f):
        return to_hat(f) * (1j * KY)
    def ddz(f):
        return to_hat(f) * (1j * KZ)
    fx = ddx(txx) + ddy(txy) + ddz(txz)
    fy = ddx(txy) + ddy(tyy) + ddz(tyz)
    fz = ddx(txz) + ddy(tyz) + ddz(tzz)
    fx, fy, fz = dealias(fx, fy, fz, mask)
    # model dissipation density proxy: <u · F>
    up, vp, wp = u, v, w
    fxp, fyp, fzp = to_phys(fx), to_phys(fy), to_phys(fz)
    eps_model = -np.mean(up * fxp + vp * fyp + wp * fzp)
    return fx, fy, fz, float(eps_model)


def force_tzt(uh, vh, wh, KX, KY, KZ, K2, n, mask, c_o, c_x, omega_c, beta, eps):
    """Sealed Phase-3 operators.
    F_diss = -curl(nu_TZT * omega),  nu_TZT = C_O * ell^2 * |omega| * Phi
    F_axis = c_X * Omega_c * P_L[ Phi * (u x omega_hat_eps) ]
    """
    ell = 2.0 * np.pi / n
    cx, cy, cz = curl_hat(uh, vh, wh, KX, KY, KZ)
    wx, wy, wz = to_phys(cx), to_phys(cy), to_phys(cz)
    wmag = np.sqrt(wx * wx + wy * wy + wz * wz)
    xi = wmag / omega_c
    phi = phi_switch(xi, beta)
    # selective viscosity
    nu_t = c_o * ell ** 2 * wmag * phi
    vx, vy, vz = nu_t * wx, nu_t * wy, nu_t * wz
    vx_h, vy_h, vz_h = to_hat(vx), to_hat(vy), to_hat(vz)
    # F_diss = - i k × Vhat
    dx, dy, dz = curl_hat(vx_h, vy_h, vz_h, KX, KY, KZ)
    fx_d, fy_d, fz_d = -dx, -dy, -dz
    # axial gate
    u, v, w = to_phys(uh), to_phys(vh), to_phys(wh)
    inv = 1.0 / np.sqrt(wmag * wmag + eps * eps)
    ox, oy, oz = wx * inv, wy * inv, wz * inv
    # u × omega_hat
    gx = v * oz - w * oy
    gy = w * ox - u * oz
    gz = u * oy - v * ox
    gx, gy, gz = phi * gx, phi * gy, phi * gz
    gx_h, gy_h, gz_h = to_hat(gx), to_hat(gy), to_hat(gz)
    gx_h, gy_h, gz_h = project(gx_h, gy_h, gz_h, KX, KY, KZ, K2)
    fx = fx_d + c_x * omega_c * gx_h
    fy = fy_d + c_x * omega_c * gy_h
    fz = fz_d + c_x * omega_c * gz_h
    fx, fy, fz = dealias(fx, fy, fz, mask)
    fxp, fyp, fzp = to_phys(fx), to_phys(fy), to_phys(fz)
    eps_model = -np.mean(u * fxp + v * fyp + w * fzp)
    return fx, fy, fz, float(eps_model)


# ---------------------------------------------------------------------------
# RHS and RK4
# ---------------------------------------------------------------------------

def rhs(uh, vh, wh, KX, KY, KZ, K2, mask, nu, n, model, params):
    uh, vh, wh = dealias(uh, vh, wh, mask)
    # convective term in rotational form: u × omega, then project
    cx, cy, cz = curl_hat(uh, vh, wh, KX, KY, KZ)
    u, v, w = to_phys(uh), to_phys(vh), to_phys(wh)
    wx, wy, wz = to_phys(cx), to_phys(cy), to_phys(cz)
    # (u × ω)
    lx = v * wz - w * wy
    ly = w * wx - u * wz
    lz = u * wy - v * wx
    lx_h, ly_h, lz_h = to_hat(lx), to_hat(ly), to_hat(lz)
    # rotational form: -(ω × u) = u × ω, projected. Pressure absorbs |u|^2/2.
    adv_x, adv_y, adv_z = project(lx_h, ly_h, lz_h, KX, KY, KZ, K2)
    # viscous
    visc_x = -nu * K2 * uh
    visc_y = -nu * K2 * vh
    visc_z = -nu * K2 * wh
    eps_nu = nu * np.mean(wx * wx + wy * wy + wz * wz)  # <ν |ω|^2> = <ν |∇u|^2> for div-free
    # model
    eps_model = 0.0
    if model == "baseline":
        mx = my = mz = 0.0
    elif model == "smagorinsky":
        mx, my, mz, eps_model = force_smagorinsky(
            uh, vh, wh, KX, KY, KZ, n, params["cs"], mask
        )
        mx, my, mz = project(mx, my, mz, KX, KY, KZ, K2)
    elif model == "tzt":
        mx, my, mz, eps_model = force_tzt(
            uh, vh, wh, KX, KY, KZ, K2, n, mask,
            params["c_o"], params["c_x"], params["omega_c"], params["beta"], params["eps"],
        )
        # axis already projected; diss is a curl so divergence-free. Project once for safety.
        mx, my, mz = project(mx, my, mz, KX, KY, KZ, K2)
    else:
        raise ValueError(model)
    return (
        adv_x + visc_x + mx,
        adv_y + visc_y + my,
        adv_z + visc_z + mz,
        float(eps_nu),
        float(eps_model),
    )


def rk4_step(uh, vh, wh, dt, ops, nu, n, model, params):
    KX, KY, KZ, K2, mask = ops
    k1x, k1y, k1z, _, _ = rhs(uh, vh, wh, KX, KY, KZ, K2, mask, nu, n, model, params)
    k2x, k2y, k2z, _, _ = rhs(
        uh + 0.5 * dt * k1x, vh + 0.5 * dt * k1y, wh + 0.5 * dt * k1z,
        KX, KY, KZ, K2, mask, nu, n, model, params,
    )
    k3x, k3y, k3z, _, _ = rhs(
        uh + 0.5 * dt * k2x, vh + 0.5 * dt * k2y, wh + 0.5 * dt * k2z,
        KX, KY, KZ, K2, mask, nu, n, model, params,
    )
    k4x, k4y, k4z, eps_nu, eps_model = rhs(
        uh + dt * k3x, vh + dt * k3y, wh + dt * k3z,
        KX, KY, KZ, K2, mask, nu, n, model, params,
    )
    uh = uh + (dt / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
    vh = vh + (dt / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y)
    wh = wh + (dt / 6.0) * (k1z + 2 * k2z + 2 * k3z + k4z)
    uh, vh, wh = project(uh, vh, wh, KX, KY, KZ, K2)
    uh, vh, wh = dealias(uh, vh, wh, mask)
    return uh, vh, wh, eps_nu, eps_model


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_case(name, model, n, t_end, dt, nu, params, out_dir, sample_every):
    X, Y, Z, KX, KY, KZ, K2, K, mask = make_grid(n)
    ops = (KX, KY, KZ, K2, mask)
    u, v, w = taylor_green(X, Y, Z)
    uh, vh, wh = to_hat(u), to_hat(v), to_hat(w)
    uh, vh, wh = project(uh, vh, wh, KX, KY, KZ, K2)
    uh, vh, wh = dealias(uh, vh, wh, mask)

    times, energies, eps_nu, eps_model, winf = [], [], [], [], []
    spectrum_at_peak = None
    peak_t = None
    t = 0.0
    step = 0
    # sample window for the reported peak spectrum
    spec_lo, spec_hi = 8.0, 9.0
    print(f"[{name}] N={n} dt={dt} t_end={t_end}")
    while t < t_end - 1e-12:
        if step % sample_every == 0:
            e = energy(uh, vh, wh, n)
            # dissipation split from a diagnostic RHS
            _, _, _, enu, emod = rhs(uh, vh, wh, KX, KY, KZ, K2, mask, nu, n, model, params)
            wi = max_vorticity(uh, vh, wh, KX, KY, KZ)
            times.append(t)
            energies.append(e)
            eps_nu.append(enu)
            eps_model.append(emod)
            winf.append(wi)
            if spec_lo <= t <= spec_hi and spectrum_at_peak is None:
                spectrum_at_peak = shell_spectrum(uh, vh, wh, K, n)
                peak_t = t
            print(f"  t={t:6.3f}  E={e:.6e}  eps_nu={enu:.4e}  eps_m={emod:.4e}  |w|_inf={wi:.3f}")
        uh, vh, wh, _, _ = rk4_step(uh, vh, wh, dt, ops, nu, n, model, params)
        t += dt
        step += 1
        if not np.isfinite(energy(uh, vh, wh, n)):
            print(f"[{name}] non-finite state at t={t:.3f}; aborting case")
            break
    # final sample
    e = energy(uh, vh, wh, n)
    _, _, _, enu, emod = rhs(uh, vh, wh, KX, KY, KZ, K2, mask, nu, n, model, params)
    wi = max_vorticity(uh, vh, wh, KX, KY, KZ)
    times.append(t_end if t >= t_end else t)
    energies.append(e)
    eps_nu.append(enu)
    eps_model.append(emod)
    winf.append(wi)
    if spectrum_at_peak is None:
        spectrum_at_peak = shell_spectrum(uh, vh, wh, K, n)
        peak_t = times[-1]
    return {
        "name": name,
        "t": np.array(times),
        "E": np.array(energies),
        "eps_nu": np.array(eps_nu),
        "eps_model": np.array(eps_model),
        "winf": np.array(winf),
        "spectrum": spectrum_at_peak,
        "peak_t": peak_t,
    }


def save_csv(results, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "t", "E", "eps_nu", "eps_model", "eps_total", "winf"])
        for r in results:
            for i in range(len(r["t"])):
                et = r["eps_nu"][i] + r["eps_model"][i]
                w.writerow([r["name"], r["t"][i], r["E"][i], r["eps_nu"][i], r["eps_model"][i], et, r["winf"][i]])


def save_plots(results, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"baseline": "#1f77b4", "smagorinsky": "#d62728", "tzt": "#2ca02c"}
    labels = {
        "baseline": "Baseline DNS",
        "smagorinsky": "Smagorinsky Cs=0.17",
        "tzt": "NS_TZT",
    }

    # dissipation rate
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for r in results:
        et = r["eps_nu"] + r["eps_model"]
        ax.plot(r["t"], et, color=colors[r["name"]], lw=2.0, label=labels[r["name"]] + " total")
        if r["name"] != "baseline":
            ax.plot(r["t"], r["eps_model"], color=colors[r["name"]], lw=1.2, ls="--", label=labels[r["name"]] + " model")
    ax.set_xlabel("t")
    ax.set_ylabel(r"dissipation rate $\varepsilon(t)$")
    ax.set_title("TGV Re=1600: dissipation rate")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "tgv_dissipation_rate.png", dpi=140)
    plt.close(fig)

    # spectrum
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for r in results:
        ek = r["spectrum"]
        k = np.arange(len(ek))
        m = (k >= 1) & (ek > 0)
        ax.loglog(k[m], ek[m], color=colors[r["name"]], lw=2.0, label=f"{labels[r['name']]} (t={r['peak_t']:.2f})")
    # Kolmogorov guide
    kk = np.array([2.0, 16.0])
    ax.loglog(kk, 0.05 * kk ** (-5.0 / 3.0), color="0.5", ls=":", label=r"$k^{-5/3}$ guide")
    ax.set_xlabel("k")
    ax.set_ylabel("E(k)")
    ax.set_title("TGV Re=1600: shell spectrum near peak")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "tgv_energy_spectrum.png", dpi=140)
    plt.close(fig)

    # max vorticity (extra, useful)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for r in results:
        ax.plot(r["t"], r["winf"], color=colors[r["name"]], lw=2.0, label=labels[r["name"]])
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\|\omega\|_\infty$")
    ax.set_title("TGV Re=1600: maximum vorticity")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "tgv_max_vorticity.png", dpi=140)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=64, help="grid (64 or 128)")
    p.add_argument("--t-end", type=float, default=10.0)
    p.add_argument("--dt", type=float, default=None, help="RK4 step; default scales with N")
    p.add_argument("--sample-every", type=int, default=None)
    p.add_argument("--out", type=str, default=".")
    args = p.parse_args()

    n = args.n
    t_end = args.t_end
    # CFL: U~1, kmax~n/3 after dealias. dt ~ 0.4 / kmax is safe for RK4 rotational form.
    dt = args.dt if args.dt is not None else min(0.01, 0.45 / (n / 3.0))
    sample_every = args.sample_every if args.sample_every is not None else max(1, int(round(0.25 / dt)))
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    nu = 1.0 / 1600.0
    common = dict(n=n, t_end=t_end, dt=dt, nu=nu, out_dir=out_dir, sample_every=sample_every)

    cases = [
        ("baseline", "baseline", {}),
        ("smagorinsky", "smagorinsky", {"cs": 0.17}),
        ("tzt", "tzt", {"c_o": 0.8, "c_x": 0.1, "omega_c": 2.0, "beta": 6.0, "eps": 1e-6}),
    ]
    results = []
    for name, model, params in cases:
        results.append(run_case(name, model, params=params, **common))

    save_csv(results, out_dir / "tgv_metrics.csv")
    save_plots(results, out_dir)
    print(f"wrote {out_dir / 'tgv_metrics.csv'}")
    print(f"wrote {out_dir / 'tgv_dissipation_rate.png'}")
    print(f"wrote {out_dir / 'tgv_energy_spectrum.png'}")


if __name__ == "__main__":
    main()
