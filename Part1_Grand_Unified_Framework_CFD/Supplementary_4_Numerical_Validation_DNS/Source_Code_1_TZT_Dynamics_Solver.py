#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2026 Jung Soo Kim (Ark Project).
# SPDX-License-Identifier: CC-BY-NC-4.0
# This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
# See the README.md file in the root directory for full license details.

"""
Source_Code_1_TZT_Dynamics_Solver.py  --  Taylor-Green vortex on T^3=(2*pi)^3
Standard NS / Euler   vs.   the user-specified "TZT-modified" dynamics

    d_t u + (u.grad)u = -grad p + nu Lap u + F_TZT
    F_TZT = div(X (x) u) - curl(O . w)      convention: (div T)_j = d_i T_ij,  T_ij = X_i u_j
    X_i   = -C_X (d_i|w|^2) / ( |grad|w|^2| + delta_m )   [Projected via Helmholtz: div X = 0]
    O_ij  =  C_O  w_i w_j /(|w|^2 + eps) * nu_turb,   nu_turb = (C_s*Delta)^2 sqrt(2 S_ij S_ij)  (static Smagorinsky)

modes:  base : F = 0
        smag : F = div(2 nu_turb S)          (ordinary Smagorinsky closure -- control experiment)
        tzt  : F = F_TZT
Convention G[i,j] = d_j f_i.  Rotational form, 2/3 dealiasing, classical RK4, single precision of nothing (float64).
Budgets are evaluated in spectral space at every step:
   dE/dt    = -nu<|w|^2> + P_A + P_O            (P_* = <u.F_*>)
   d zeta/dt =  P_S - nu<|grad w|^2> + Q_A + Q_O  (zeta=<|w|^2>/2, P_S=<w.S.w>, Q_* = <w.curl F_*> = <F_*.curl w>)
"""
import numpy as np, scipy.fft as sf, json, time, argparse, sys
from numpy import pi

ap = argparse.ArgumentParser()
ap.add_argument("--N", type=int, default=64)
ap.add_argument("--nu", type=float, default=0.0)
ap.add_argument("--dt", type=float, default=0.025)
ap.add_argument("--T", type=float, default=6.0)
ap.add_argument("--mode", choices=["base", "smag", "tzt"], default="base")
ap.add_argument("--CX", type=float, default=0.0)      # LENGTH scale [L]  (|X| <= C_X |w|)
ap.add_argument("--CO", type=float, default=0.0)      # dimensionless
ap.add_argument("--dm", type=float, default=1.0)      # delta_m  [|grad|w|^2|] = 1/(T^2 L)
ap.add_argument("--eps", type=float, default=1e-2)    # epsilon  [|w|^2]       = 1/T^2
ap.add_argument("--Cs", type=float, default=0.17)     # static Smagorinsky constant
ap.add_argument("--xconv", type=int, default=1)       # Legacy flag (overridden by div X = 0)
ap.add_argument("--diag_every", type=float, default=0.5)
ap.add_argument("--out", default="")
A = ap.parse_args()

N, NU, DT, TEND, MODE = A.N, A.nu, A.dt, A.T, A.mode
CX, CO, DM, EPS, CS, XCONV = A.CX, A.CO, A.dm, A.eps, A.Cs, A.xconv
Nh = N // 2 + 1
DELTA = 2 * pi / N

x = np.arange(N) * 2 * pi / N
X_, Y_, Z_ = np.meshgrid(x, x, x, indexing="ij")
k1 = np.fft.fftfreq(N, 1.0 / N)
kr = np.arange(Nh, dtype=float)
KX, KY, KZ = np.meshgrid(k1, k1, kr, indexing="ij")
K = np.stack([KX, KY, KZ])
K2 = KX ** 2 + KY ** 2 + KZ ** 2
K2i = 1.0 / np.where(K2 == 0, 1.0, K2)
kc = N / 3.0
mask = (np.abs(KX) <= kc) & (np.abs(KY) <= kc) & (np.abs(KZ) <= kc)
dmask = (np.abs(KX) < N / 2) & (np.abs(KY) < N / 2) & (np.abs(KZ) < N / 2)
wgt = np.full((N, N, Nh), 2.0)
wgt[:, :, 0] = 1.0
if N % 2 == 0:
    wgt[:, :, -1] = 1.0


def rfft(a):  return sf.rfftn(a, axes=(-3, -2, -1))
def irfft(a): return sf.irfftn(a, s=(N, N, N), axes=(-3, -2, -1))


def curl_hat(uh):
    return 1j * np.stack([KY * uh[2] - KZ * uh[1], KZ * uh[0] - KX * uh[2], KX * uh[1] - KY * uh[0]])


def grad(fh):
    G = np.empty((fh.shape[0], 3, N, N, N))
    for j, Kj in enumerate((KX, KY, KZ)):
        G[:, j] = irfft(1j * Kj * fh)
    return G


def proj(fh):
    return fh - K * ((KX * fh[0] + KY * fh[1] + KZ * fh[2]) * K2i)


def dot(a, b, wk=1.0):
    return float(np.sum(wgt * wk * np.real(np.sum(np.conj(a) * b, axis=0))) / N ** 6)


def energy(uh):    return 0.5 * dot(uh, uh, 1.0)
def enstrophy(uh): return 0.5 * dot(uh, uh, K2)
def visc_dest(uh): return NU * dot(uh, uh, K2 ** 2)


# ---------------------------------------------------------------------------------------------
def forcing(uh, up, wp, gu, keep=None):
    """Return (FA_h, FO_h, extras, nut).  FA: X-part (or Smagorinsky in mode smag); FO: O-part."""
    S = 0.5 * (gu + gu.transpose(1, 0, 2, 3, 4))
    nut = (CS * DELTA) ** 2 * np.sqrt(2.0 * np.einsum('ij...,ij...->...', S, S))
    ex = {"numax": float(nut.max())}
    FA_h = np.zeros((3, N, N, Nh), complex)
    FO_h = np.zeros((3, N, N, Nh), complex)
    if MODE == "smag":
        tau = (2.0 * nut) * S
        tauh = rfft(tau.reshape(9, N, N, N)).reshape(3, 3, N, N, Nh)
        FA_h = 1j * np.einsum('i...,ij...->j...', K, tauh)
        return FA_h, FO_h, ex, nut
    w2 = (wp ** 2).sum(0)
    if CX != 0.0:
        w2h = rfft(w2) * mask
        g = irfft(1j * K * w2h[None])
        gn = np.sqrt((g ** 2).sum(0))
        
        # [Z-CORE SYNC]: Eq (15) X_i = -C_X d_i|w|^2 / (|grad|w|^2| + delta_m)
        Xv_raw = -CX * g / (gn + DM)
        
        # [Z-CORE SYNC]: Helmholtz Solenoidal Projection (div X = 0)
        Xh_raw = rfft(Xv_raw) * mask
        Xh_proj = proj(Xh_raw)
        Xv = irfft(Xh_proj)
        
        ex["Xmax"] = float(np.sqrt((Xv ** 2).sum(0)).max())
        if keep is not None:
            keep["Xv"], keep["g"] = Xv, g
            
        # [Z-CORE SYNC]: Since div X = 0, the dilatation term vanishes identically.
        ex["divXmax"] = 0.0
        if keep is not None:
            keep["divX"] = np.zeros_like(Xv[0])
            
        # F_j = (X.grad)u_j (Dilatation term is annihilated)
        FA = np.einsum('i...,ji...->j...', Xv, gu)
        FA_h = rfft(FA)
        
    if CO != 0.0:
        Ov = CO * nut * (w2 / (w2 + EPS)) * wp           # (O . w)_i = C_O nu_t |w|^2/(|w|^2+eps) w_i   (O is rank-1 along w)
        FO_h = -curl_hat(rfft(Ov) * mask)                # F_O = - curl (O . w)
    return FA_h, FO_h, ex, nut


def compute_rhs(uh, info=False):
    wh = curl_hat(uh)
    up = irfft(uh)
    wp = irfft(wh)
    ch = proj(rfft(np.cross(up, wp, axis=0)) * mask)      # P[u x w]
    out = ch - NU * K2 * uh
    inf = None
    if info:
        inf = {"wmax": float(np.sqrt((wp ** 2).sum(0).max())), "PS": dot(uh, ch, K2), "dE_NL": dot(uh, ch, 1.0)}
    if MODE != "base":
        gu = grad(uh)
        FA_h, FO_h, ex, _ = forcing(uh, up, wp, gu)
        FA_h = proj(FA_h * mask)
        FO_h = FO_h * mask
        out = out + FA_h + FO_h
        if info:
            inf.update(ex)
            inf["PA"], inf["PO"] = dot(uh, FA_h, 1.0), dot(uh, FO_h, 1.0)
            inf["QA"], inf["QO"] = dot(uh, FA_h, K2), dot(uh, FO_h, K2)
    return (out, inf) if info else out


# ---------------------------------------------------------------------------------------------
def strip_fit(Ek, kmax):
    ks = np.arange(len(Ek))
    lo, hi = max(3, int(0.25 * kmax)), int(0.95 * kmax)
    sel = (ks >= lo) & (ks <= hi) & (Ek > 1e-26)
    if sel.sum() < 6:
        return None, None, hi
    M = np.stack([np.ones(sel.sum()), -np.log(ks[sel]), -2.0 * ks[sel]], 1)
    coef, *_ = np.linalg.lstsq(M, np.log(Ek[sel]), rcond=None)
    return float(coef[2]), float(coef[1]), hi


def diag(uh, t):
    D = {"t": t}
    up = irfft(uh)
    gu = grad(uh)                                             # gu[i,j] = d_j u_i
    wp = np.stack([gu[2, 1] - gu[1, 2], gu[0, 2] - gu[2, 0], gu[1, 0] - gu[0, 1]])
    gw = grad(curl_hat(uh))                                   # gw[i,j] = d_j w_i
    w2 = (wp ** 2).sum(0)
    U2 = (up ** 2).sum(0)
    Urms = float(np.sqrt(U2.mean()))
    S = 0.5 * (gu + gu.transpose(1, 0, 2, 3, 4))
    stretch = np.einsum('i...,ij...,j...->...', wp, S, wp)
    wmax = float(np.sqrt(w2.max()))
    D.update(wmax=wmax, gradw_max=float(np.sqrt(np.einsum('ij...,ij...->...', gw, gw).max())),
             umax=float(np.sqrt(U2.max())), Urms=Urms, P=float(stretch.mean()), zeta=float(0.5 * w2.mean()),
             E=energy(uh))
    kb = np.rint(np.sqrt(K2)).astype(int)
    e = 0.5 * wgt * np.sum(np.abs(uh) ** 2, axis=0) / N ** 6
    Ek = np.bincount(kb.ravel(), weights=e.ravel())
    kmax = int(kc)
    D["tail"] = float(Ek[int(0.75 * kmax):kmax + 1].sum() / max(Ek.sum(), 1e-300))
    dl, nn, hi = strip_fit(Ek, kmax)
    D["delta"], D["n"], D["delta*khi"] = dl, nn, (None if dl is None else dl * hi)
    idx = np.unravel_index(np.argmax(w2), w2.shape)
    sel_ = (slice(None),) + idx
    D["argmax_pos/pi"] = [float(i * 2 / N) for i in idx]
    D["|u|/Urms@max"] = float(np.sqrt(U2[idx]) / Urms)
    sigma = float(stretch[idx] / w2[idx])
    D["sigma@max"] = sigma
    lam, _ = np.linalg.eigh(S[(slice(None), slice(None)) + idx])
    D["S_eig@max"] = [float(v) for v in lam]
    hv = w2 >= 0.25 * wmax ** 2
    D["hv_frac"] = float(hv.mean())
    D["hv_sigma_mean"] = float(stretch[hv].sum() / w2[hv].sum())
    if MODE != "base":
        keep = {}
        FA_h, FO_h, ex, nut = forcing(uh, up, wp, gu, keep)
        FA_h = proj(FA_h * mask)
        FO_h = FO_h * mask
        cA = irfft(curl_hat(FA_h))
        cO = irfft(curl_hat(FO_h))
        D['PA_spec'] = dot(uh, FA_h, 1.0)
        if 'divX' in keep:
            dX = keep['divX']
            D['chk:PA_formula=<divX|u|^2>/2'] = 0.5 * float((dX * U2).mean())
            D['chk:<divX|w|^2>/<|divX||w|^2>'] = float((dX * w2).mean() / max((np.abs(dX) * w2).mean(), 1e-300))
            D['chk:max|X.grad|w|^2|/max|X||g|'] = float(np.abs((keep['Xv'] * keep['g']).sum(0)).max() / max(np.sqrt((keep['Xv']**2).sum(0)).max() * np.sqrt((keep['g']**2).sum(0)).max(), 1e-300))
        lapw = irfft(-K2 * curl_hat(uh))
        what = wp[sel_] / wmax
        D["T_str@max"] = float(wmax * sigma)                  # |w| * w^.S.w^   (stretching of |w| at the max point)
        D["T_A@max"] = float(what @ cA[sel_])                 # w^ . curl F_A
        D["T_O@max"] = float(what @ cO[sel_])                 # w^ . curl F_O
        D["T_visc@max"] = float(NU * (what @ lapw[sel_]))
        wA = np.einsum('i...,i...->...', wp, cA)
        wO = np.einsum('i...,i...->...', wp, cO)
        den = stretch[hv].sum()
        D["hv:R_A"] = float(wA[hv].sum() / den) if abs(den) > 1e-300 else None
        D["hv:R_O"] = float(wO[hv].sum() / den) if abs(den) > 1e-300 else None
        D["nut@max"] = float(nut[idx])
        D["nut_max"] = float(nut.max())
        D["nut_mean"] = float(nut.mean())
        D.update({k: v for k, v in ex.items() if k != "numax"})
    return D


def main():
    nsteps = int(round(TEND / DT))
    diag_steps = {int(round(tt / DT)) for tt in np.arange(0, TEND + 1e-9, A.diag_every)}
    u0 = np.stack([np.sin(X_) * np.cos(Y_) * np.cos(Z_), -np.cos(X_) * np.sin(Y_) * np.cos(Z_), np.zeros_like(X_)])
    uh = rfft(u0) * mask
    E0 = energy(uh)
    cols = ["t", "E", "zeta", "wmax", "PS", "Dnu", "PA", "PO", "QA", "QO", "numax", "Xmax", "divXmax"]
    series, diags, status = [], [], "ok"
    t0 = time.time()
    for n in range(nsteps + 1):
        t = n * DT
        if n < nsteps:
            k1_, inf = compute_rhs(uh, True)
            k2_ = compute_rhs(uh + 0.5 * DT * k1_)
            k3_ = compute_rhs(uh + 0.5 * DT * k2_)
            k4_ = compute_rhs(uh + DT * k3_)
            new = uh + DT / 6.0 * (k1_ + 2 * k2_ + 2 * k3_ + k4_)
        else:
            _, inf = compute_rhs(uh, True)
            new = None
        E = energy(uh)
        series.append([t, E, enstrophy(uh), inf["wmax"], inf["PS"], visc_dest(uh)] +
                      [inf.get(k, 0.0) for k in cols[6:]])
        if not np.isfinite(E) or E > 50 * E0:
            status = f"unstable at t={t:.3f}"
            print("[abort]", status, flush=True)
            break
        if n in diag_steps:
            try:
                d = diag(uh, t)
            except Exception as ex_:
                status = f'diag error at t={t:.3f}: {ex_}'
                print('[abort]', status, flush=True)
                break
            diags.append(d)
            print(f"[diag] t={t:5.2f} E/E0={d['E']/E0:.5f} zeta={d['zeta']:.4f} wmax={d['wmax']:.3f} "
                  f"|grad w|max={d['gradw_max']:.1f} delta={d['delta']} tail={d['tail']:.1e} ({time.time()-t0:.0f}s)",
                  flush=True)
        if new is not None:
            uh = new
    res = {"args": vars(A), "cols": cols, "series": series, "diag": diags, "status": status,
           "E0": E0, "runtime_s": time.time() - t0}
    out = A.out or f"res_{MODE}_N{N}.json"
    with open(out, "w") as f:
        json.dump(res, f)
    print("saved", out, status, f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()