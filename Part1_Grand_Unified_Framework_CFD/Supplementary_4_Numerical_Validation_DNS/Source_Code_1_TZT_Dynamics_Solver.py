#!/usr/bin/env python3
"""
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

Source_Code_1_TZT_Dynamics_Solver.py  --  v3 (localized TZT activation, "Fractal Valve")

Pseudo-spectral solver for the Taylor-Green vortex on T^3 = (0, 2 pi)^3: float64, Fourier
derivatives, 2/3 dealiasing, Leray projection, classical RK4.

Modes (--mode)
  base : incompressible Navier-Stokes / Euler, no extra forcing
  smag : static Smagorinsky closure (control experiment), applied everywhere
  tzt  : TZT forcing  F = div(X (x) u) - curl(O . omega)
           X   = P[ -C_X grad|omega|^2 / (|grad|omega|^2| + delta_m) ]   (solenoidal, so
                 div(X (x) u) = (X . grad) u)
           O.omega = C_O nu_t |omega|^2/(|omega|^2 + eps) omega,
           nu_t = (C_s Delta)^2 sqrt(2 S:S),  Delta = 2 pi / N

v3 -- synchronized with Lean 4 v3 (Supplementary_2: Theorem0, 1, 5, 6)
  The TZT forcing is no longer applied on the whole domain.  It is multiplied by the
  Fractal Valve chi(x) in [0, 1], supported in the critical region K:
      --valve omega  : K = {x : |omega(x)|^2 > Lambda^2}      (= Lean `CriticalRegion Lambda`)
      --valve strain : K = {x : |S(x)|^2 > Lambda_S^2},  |S| = sqrt(2 S_ij S_ij)
      --valve either : union of the two
      --valve global : K = T^3 (v2 behaviour, for comparison only)
  chi is the sharp indicator of K (--valve_width 0, default), or a C^1 ramp that is 0 for
  |q| <= Lambda and 1 for |q| >= (1 + w) Lambda, so supp chi is inside K in every case.
  Outside K the forcing is exactly zero, so the momentum equation there is classical
  Navier-Stokes / Euler; the pressure, as in any incompressible flow, stays global.
  Because chi F is no longer solenoidal, both forcing pieces are Leray-projected.
  Thresholds: Lambda = --Lambda if given (>= 0), else --Lambda_rel * max|omega(t=0)|;
              Lambda_S = --LambdaS if given (>= 0), else --Lambda_rel * max|S(t=0)|.

  Reading the diagnostics against Lean.  Lean `vorticity_ceiling` says: a state that
  satisfies the localized postulate (S(v) = 0 on a core containing CriticalRegion Lambda)
  has |omega| <= Lambda.  This solver does not impose that postulate as a constraint; it
  applies the TZT forcing on K.  Whether a simulated flow stays below Lambda is therefore an
  output of the run (wmax/Lambda, crit_frac), and core_res measures how far the flow is from
  the postulate S(v) = 0 on the active region.

Other v3 changes: divXmax is measured (v2 stored 0.0); the series gains the columns
valve_frac and core_res; the JSON gains version, Lambda, LambdaS, wmax0, Smax0.
"""

import argparse
import json
import time

import numpy as np
import scipy.fft as sf
from numpy import pi

ap = argparse.ArgumentParser(description="Taylor-Green dynamics: NS/Euler vs TZT (v3, localized)")
ap.add_argument("--N", type=int, default=64, help="grid points per direction")
ap.add_argument("--nu", type=float, default=0.0, help="kinematic viscosity")
ap.add_argument("--dt", type=float, default=0.025, help="time step")
ap.add_argument("--T", type=float, default=6.0, help="end time")
ap.add_argument("--mode", choices=["base", "smag", "tzt"], default="base")
ap.add_argument("--CX", type=float, default=0.0, help="C_X (length scale of X)")
ap.add_argument("--CO", type=float, default=0.0, help="C_O (dimensionless)")
ap.add_argument("--dm", type=float, default=1.0, help="delta_m regulariser in X")
ap.add_argument("--eps", type=float, default=1e-2, help="epsilon regulariser in O")
ap.add_argument("--Cs", type=float, default=0.17, help="Smagorinsky constant C_s")
ap.add_argument("--xconv", type=int, default=1, help="legacy flag (unused)")
ap.add_argument("--diag_every", type=float, default=0.5, help="diagnostic interval in time")
ap.add_argument("--out", type=str, default="", help="output path (default res_{mode}_N{N}.json)")
# v3: Fractal Valve
ap.add_argument("--valve", choices=["omega", "strain", "either", "global"], default="omega",
                help="where the TZT forcing is active (default: Lean CriticalRegion on |omega|)")
ap.add_argument("--Lambda", type=float, default=-1.0,
                help="absolute |omega| threshold Lambda (< 0: use --Lambda_rel)")
ap.add_argument("--LambdaS", type=float, default=-1.0,
                help="absolute |S| threshold Lambda_S (< 0: use --Lambda_rel)")
ap.add_argument("--Lambda_rel", type=float, default=1.0,
                help="Lambda = Lambda_rel*max|omega(0)|, Lambda_S = Lambda_rel*max|S(0)|")
ap.add_argument("--valve_width", type=float, default=0.0,
                help="0: sharp indicator of K; w > 0: C^1 ramp from Lambda to (1+w)Lambda")
args = ap.parse_args()

# ----------------------------------------------------------------------------- grid
N = args.N
Nh = N // 2 + 1
DELTA = 2 * pi / N
xg = np.arange(N) * 2 * pi / N
GX, GY, GZ = np.meshgrid(xg, xg, xg, indexing="ij")
k1 = sf.fftfreq(N, 1.0 / N)
kr = np.arange(Nh)
KX, KY, KZ = np.meshgrid(k1, k1, kr, indexing="ij")
K = np.array([KX, KY, KZ])
K2 = KX**2 + KY**2 + KZ**2
K2i = 1.0 / np.where(K2 == 0, 1.0, K2)
kc = N / 3.0
mask = ((np.abs(KX) <= kc) & (np.abs(KY) <= kc) & (np.abs(KZ) <= kc)).astype(float)
wgt = np.full(KZ.shape, 2.0)            # Parseval weights for the half spectrum
wgt[..., 0] = 1.0
if N % 2 == 0:
    wgt[..., N // 2] = 1.0
KB = np.rint(np.sqrt(K2)).astype(int)   # shell index
KMAX = N // 3

LAM = 0.0                               # set from the initial condition below
LAMS = 0.0


# ----------------------------------------------------------------------------- spectral helpers
def rfft(a):
    return sf.rfftn(a, axes=(-3, -2, -1))


def irfft(a):
    return sf.irfftn(a, s=(N, N, N), axes=(-3, -2, -1))


def curl_hat(uh):
    return np.array([1j * (KY * uh[2] - KZ * uh[1]),
                     1j * (KZ * uh[0] - KX * uh[2]),
                     1j * (KX * uh[1] - KY * uh[0])])


def grad(fh):
    """fh: (m, N, N, Nh) -> (m, 3, N, N, N), G[a, j] = d_j f_a (no dealiasing)."""
    return np.array([[irfft(1j * K[j] * fh[a]) for j in range(3)] for a in range(fh.shape[0])])


def proj(fh):
    """Leray projection P f = f - k (k . f)/|k|^2."""
    kdotf = KX * fh[0] + KY * fh[1] + KZ * fh[2]
    return fh - K * (kdotf * K2i)


def dot(a, b, wk=1.0):
    """Grid-averaged inner product <a, b> from half-spectrum coefficients."""
    return float(np.sum(wgt * wk * np.real(np.sum(np.conj(a) * b, axis=0)))) / N**6


def energy(uh):
    return 0.5 * dot(uh, uh, 1.0)


def enstrophy(uh):
    return 0.5 * dot(uh, uh, K2)


def visc_dest(uh):
    """nu <|grad omega|^2> (enstrophy-budget viscous term), stored as 'Dnu'."""
    return args.nu * dot(uh, uh, K2**2)


def strain_of(gu):
    return 0.5 * (gu + gu.transpose(1, 0, 2, 3, 4))


# ----------------------------------------------------------------------------- v3: Fractal Valve
def gate(q2, L, width):
    """Valve on a squared magnitude q2 with threshold L >= 0.
    width = 0 : sharp indicator of {q2 > L^2}   (Lean CriticalRegion: Lambda^2 < |omega|^2)
    width > 0 : C^1 smoothstep in |q| from 0 at L to 1 at (1+width) L; identically 0 on
                {|q| <= L}, so the support stays inside {q2 > L^2}."""
    if width <= 0.0 or L <= 0.0:
        return (q2 > L * L).astype(float)
    s = np.clip((np.sqrt(q2) - L) / (width * L), 0.0, 1.0)
    return s * s * (3.0 - 2.0 * s)


def valve(w2, S2):
    """chi(x) in [0, 1]; S2 = 2 S_ij S_ij."""
    if args.valve == "global":
        return np.ones_like(w2)
    chi = np.zeros_like(w2)
    if args.valve in ("omega", "either"):
        chi = np.maximum(chi, gate(w2, LAM, args.valve_width))
    if args.valve in ("strain", "either"):
        chi = np.maximum(chi, gate(S2, LAMS, args.valve_width))
    return chi


def core_residual(chi, S2):
    """sqrt(<2 S:S>_chi) / Lambda: distance from the Lean core condition S(v) = 0 on K."""
    m = float(np.sum(chi))
    if m <= 0.0 or LAM <= 0.0:
        return 0.0
    return float(np.sqrt(np.sum(chi * S2) / m) / LAM)


# ----------------------------------------------------------------------------- forcing
def forcing(uh, up, wp, gu, keep=None):
    """Returns (FA_h, FO_h, ex, nut); FA_h, FO_h are spectral (unmasked, unprojected)."""
    S = strain_of(gu)
    S2 = 2.0 * np.sum(S * S, axis=(0, 1))
    nut = (args.Cs * DELTA) ** 2 * np.sqrt(S2)
    ex = {"numax": float(nut.max())}
    if args.mode == "smag":
        tauh = rfft(2.0 * nut[None, None] * S)
        F = np.array([1j * sum(K[i] * tauh[i, j] for i in range(3)) for j in range(3)])
        return F, np.zeros_like(F), ex, nut

    w2 = np.sum(wp * wp, axis=0)
    FA = np.zeros((3, N, N, N))
    if args.CX != 0.0:
        g = grad((rfft(w2) * mask)[None])[0]                 # g_i = d_i |omega|^2
        gn = np.sqrt(np.sum(g * g, axis=0))
        Xraw = -args.CX * g / (gn + args.dm)
        Xh = proj(rfft(Xraw) * mask)
        X = irfft(Xh)
        divX = irfft(1j * (KX * Xh[0] + KY * Xh[1] + KZ * Xh[2]))
        ex["Xmax"] = float(np.sqrt(np.sum(X * X, axis=0)).max())
        ex["divXmax"] = float(np.abs(divX).max())             # v3: measured
        FA = np.array([sum(X[i] * gu[j, i] for i in range(3)) for j in range(3)])   # (X.grad)u
        if keep is not None:
            keep.update(X=X, g=g, divX=divX)
    FO = np.zeros((3, N, N, N))
    if args.CO != 0.0:
        Ov = args.CO * nut * w2 / (w2 + args.eps) * wp        # O . omega
        FO = irfft(-curl_hat(rfft(Ov) * mask))

    # v3: the TZT forcing acts only on the critical region K
    chi = valve(w2, S2)
    FA = chi[None] * FA
    FO = chi[None] * FO
    ex["valve_frac"] = float(np.mean(chi > 0.0))
    ex["chi_mean"] = float(chi.mean())
    ex["core_res"] = core_residual(chi, S2)
    if keep is not None:
        keep.update(chi=chi, S2=S2)
    return rfft(FA), rfft(FO), ex, nut


# ----------------------------------------------------------------------------- right-hand side
def compute_rhs(uh, info=False):
    wh = curl_hat(uh)
    up = irfft(uh)
    wp = irfft(wh)
    uxw = np.array([up[1] * wp[2] - up[2] * wp[1],
                    up[2] * wp[0] - up[0] * wp[2],
                    up[0] * wp[1] - up[1] * wp[0]])
    ch = proj(rfft(uxw) * mask)                               # P[u x omega], rotational form
    rhs = ch - args.nu * K2 * uh
    if args.mode != "base":
        FAh, FOh, ex, nut = forcing(uh, up, wp, grad(uh))
        FAh = proj(FAh * mask)
        FOh = proj(FOh * mask)                                # v3: chi F_O is not solenoidal
        rhs = rhs + FAh + FOh
    if not info:
        return rhs, None
    inf = {"wmax": float(np.sqrt(np.sum(wp * wp, axis=0).max())),
           "PS": dot(uh, ch, K2),
           "dE_NL": dot(uh, ch, 1.0)}
    if args.mode != "base":
        inf.update(ex)
        inf.update(PA=dot(uh, FAh, 1.0), PO=dot(uh, FOh, 1.0),
                   QA=dot(uh, FAh, K2), QO=dot(uh, FOh, K2))
    return rhs, inf


# ----------------------------------------------------------------------------- diagnostics
def strip_fit(Ek, kmax):
    """Least-squares fit  log E_k = a - n log k - 2 delta k  on lo <= k <= hi."""
    lo = max(3, int(np.floor(0.25 * kmax)))
    hi = int(np.floor(0.95 * kmax))
    ks = np.arange(len(Ek))
    sel = (ks >= lo) & (ks <= hi) & (Ek > 1e-26)
    if sel.sum() < 6:
        return None, None, hi
    k = ks[sel].astype(float)
    A = np.column_stack([np.ones_like(k), -np.log(k), -2.0 * k])
    coef = np.linalg.lstsq(A, np.log(Ek[sel]), rcond=None)[0]
    return float(coef[2]), float(coef[1]), hi


def diag(uh, t):
    D = {"t": float(t)}
    wh = curl_hat(uh)
    up = irfft(uh)
    wp = irfft(wh)
    gu = grad(uh)
    gw = grad(wh)
    S = strain_of(gu)
    w2 = np.sum(wp * wp, axis=0)
    U2 = np.sum(up * up, axis=0)
    stretch = np.einsum("i...,ij...,j...->...", wp, S, wp)   # omega . S omega
    wmax = float(np.sqrt(w2.max()))
    Urms = float(np.sqrt(U2.mean()))
    D.update(wmax=wmax,
             gradw_max=float(np.sqrt(np.sum(gw * gw, axis=(0, 1))).max()),
             umax=float(np.sqrt(U2.max())), Urms=Urms,
             P=float(stretch.mean()), zeta=float(0.5 * w2.mean()), E=energy(uh))
    # spectrum, tail, exponential-cutoff fit
    e3 = 0.5 * wgt * np.sum(np.abs(uh) ** 2, axis=0) / N**6
    Ek = np.bincount(KB.ravel(), weights=e3.ravel())
    lo_t = int(np.floor(0.75 * KMAX))
    D["tail"] = float(Ek[lo_t:KMAX + 1].sum() / max(Ek.sum(), 1e-300))
    delta, n_exp, hi = strip_fit(Ek, KMAX)
    D["delta"], D["n"] = delta, n_exp
    D["delta*khi"] = None if delta is None else delta * hi
    # peak-vorticity point
    idx = np.unravel_index(int(np.argmax(w2)), w2.shape)
    sl = (slice(None),) + idx
    sigma = float(stretch[idx] / w2[idx])
    D["argmax_pos/pi"] = [2.0 * i / N for i in idx]
    D["|u|/Urms@max"] = float(np.sqrt(U2[idx]) / Urms)
    D["sigma@max"] = sigma
    D["S_eig@max"] = [float(e) for e in np.linalg.eigvalsh(S[(slice(None), slice(None)) + idx])]
    # high-vorticity region
    hv = w2 >= 0.25 * wmax**2
    D["hv_frac"] = float(hv.mean())
    D["hv_sigma_mean"] = float(stretch[hv].sum() / w2[hv].sum())
    # v3: Lean CriticalRegion diagnostics (all modes)
    D["Lambda"], D["LambdaS"] = LAM, LAMS
    D["crit_frac"] = float(np.mean(w2 > LAM * LAM))
    D["wmax/Lambda"] = wmax / LAM if LAM > 0 else None

    if args.mode != "base":
        keep = {}
        FAh, FOh, ex, nut = forcing(uh, up, wp, gu, keep)
        FA = proj(FAh * mask)
        FO = proj(FOh * mask)
        cA = irfft(curl_hat(FA))
        cO = irfft(curl_hat(FO))
        D["PA_spec"] = dot(uh, FA, 1.0)
        if "divX" in keep:
            divX, X, g = keep["divX"], keep["X"], keep["g"]
            D["chk:PA_formula=<divX|u|^2>/2"] = float(0.5 * np.mean(divX * U2))
            D["chk:<divX|w|^2>/<|divX||w|^2>"] = float(
                np.mean(divX * w2) / max(np.mean(np.abs(divX) * w2), 1e-300))
            Xn = np.sqrt(np.sum(X * X, axis=0))
            gn = np.sqrt(np.sum(g * g, axis=0))
            D["chk:max|X.grad|w|^2|/max|X||g|"] = float(
                np.abs(np.sum(X * g, axis=0)).max() / max(Xn.max() * gn.max(), 1e-300))
        what = wp[sl] / wmax
        lapw = irfft(-K2 * wh)
        D["T_str@max"] = float(np.sqrt(w2[idx]) * sigma)
        D["T_A@max"] = float(what @ cA[sl])
        D["T_O@max"] = float(what @ cO[sl])
        D["T_visc@max"] = float(args.nu * (what @ lapw[sl]))
        den = float(stretch[hv].sum())
        D["hv:R_A"] = float(np.sum(np.sum(wp * cA, axis=0)[hv]) / den) if abs(den) > 1e-300 else None
        D["hv:R_O"] = float(np.sum(np.sum(wp * cO, axis=0)[hv]) / den) if abs(den) > 1e-300 else None
        D["nut@max"] = float(nut[idx])
        D["nut_max"] = float(nut.max())
        D["nut_mean"] = float(nut.mean())
        if "Xmax" in ex:
            D["Xmax"], D["divXmax"] = ex["Xmax"], ex["divXmax"]
        if args.mode == "tzt":
            chi = keep["chi"]
            D["valve"] = args.valve
            D["valve_frac"] = ex["valve_frac"]
            D["chi_mean"] = ex["chi_mean"]
            D["valve@max"] = float(chi[idx])
            D["core_res"] = ex["core_res"]
            D["P_in"] = float(np.mean(chi * stretch))         # stretching inside the valve
            D["P_out"] = float(np.mean((1.0 - chi) * stretch))  # classical region
    return D


# ----------------------------------------------------------------------------- run
u0 = np.array([np.sin(GX) * np.cos(GY) * np.cos(GZ),
               -np.cos(GX) * np.sin(GY) * np.cos(GZ),
               np.zeros_like(GX)])
uh = mask * rfft(u0)
E0 = energy(uh)

_w0 = irfft(curl_hat(uh))
_S0 = strain_of(grad(uh))
wmax0 = float(np.sqrt(np.sum(_w0 * _w0, axis=0).max()))
Smax0 = float(np.sqrt(2.0 * np.sum(_S0 * _S0, axis=(0, 1))).max())
LAM = args.Lambda if args.Lambda >= 0 else args.Lambda_rel * wmax0
LAMS = args.LambdaS if args.LambdaS >= 0 else args.Lambda_rel * Smax0
if args.mode == "tzt":
    print(f"[valve] mode={args.valve}  Lambda={LAM:.6g}  Lambda_S={LAMS:.6g}  "
          f"(max|omega0|={wmax0:.6g}, max|S0|={Smax0:.6g}, width={args.valve_width})")

nsteps = int(round(args.T / args.dt))
dsteps = {int(round(tt / args.dt)) for tt in np.arange(0.0, args.T + 1e-9, args.diag_every)}
cols = ["t", "E", "zeta", "wmax", "PS", "Dnu", "PA", "PO", "QA", "QO", "numax", "Xmax",
        "divXmax", "valve_frac", "core_res"]
series, diags = [], []
status = "ok"
t_start = time.time()

for n in range(nsteps + 1):
    t = n * args.dt
    k1, inf = compute_rhs(uh, info=True)
    if n < nsteps:
        k2, _ = compute_rhs(uh + 0.5 * args.dt * k1)
        k3, _ = compute_rhs(uh + 0.5 * args.dt * k2)
        k4, _ = compute_rhs(uh + args.dt * k3)
        uh_new = uh + (args.dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    E = energy(uh)
    series.append([t, E, enstrophy(uh), inf["wmax"], inf["PS"], visc_dest(uh)]
                  + [inf.get(c, 0.0) for c in cols[6:]])
    if not np.isfinite(E) or E > 50.0 * E0:
        status = f"unstable at t={t:.4f}"
        print(f"[abort] {status}")
        break
    if n in dsteps:
        try:
            d = diag(uh, t)
        except Exception as exc:                              # keep the run's output
            status = f"diag error at t={t:.4f}: {exc}"
            print(f"[abort] {status}")
            break
        diags.append(d)
        line = (f"t={t:6.3f}  E={d['E']:.6e}  wmax={d['wmax']:.4f}  P={d['P']:.4e}  "
                f"crit_frac={d['crit_frac']:.4f}")
        if args.mode == "tzt":
            line += f"  valve_frac={d['valve_frac']:.4f}  core_res={d['core_res']:.3f}"
        print(line)
    if n < nsteps:
        uh = uh_new

out = args.out or f"res_{args.mode}_N{N}.json"
result = {"args": vars(args), "cols": cols, "series": series, "diag": diags, "status": status,
          "E0": E0, "runtime_s": time.time() - t_start,
          "version": "v3-localized", "Lambda": LAM, "LambdaS": LAMS,
          "wmax0": wmax0, "Smax0": Smax0}
with open(out, "w") as fh:
    json.dump(result, fh, indent=1)
print(f"[done] status={status}  wrote {out}")
