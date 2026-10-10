#!/usr/bin/env python3
"""
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

CAS_Proof1_Helmholtz_Bridge.py  --  Proof 1: Helmholtz Bridge (H1), v3 (localized core)

Lean 4 counterpart: Supplementary_2_Formal_Verification_Lean4/Theorem1_Helmholtz_Bridge.lean
                    (+ Theorem0 definitions, Theorem5 witnesses)

v3 synchronisation with Lean
  * Conventions are Lean's:  (grad v)_ij = d_j v_i,  S = sym, Omega = skw,
    omega = axial(grad v) = curl v,  |a|^2 = a . a.
  * The zero-point core  sym G0 = lambda g  is imposed ONLY on the critical region
        CriticalRegion(Lambda) = { x : Lambda^2 < |omega(x)|^2 }        (strict)
    encoded here by the indicator  chi = Piecewise((1, Lambda^2 < |omega|^2), (0, True))
    and the localized condition  chi * (sym G0 - lambda g) = 0.
  * (R), Z == 0 and div v = 0 hold everywhere, as in Lean `TZTPostulates s lam K`.
  * Every check below is an executed SymPy computation guarded by an assert.
"""

from sympy import (Function, Matrix, Piecewise, Rational, cos, diff, eye, expand, pi,
                   simplify, sin, solve, symbols)

half = Rational(1, 2)
x, y, z = symbols("x y z", real=True)
coords = (x, y, z)


# --------------------------------------------------------------------------- helpers
def header(t):
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)


def is_zero(e):
    """True iff the scalar / every matrix entry simplifies to 0."""
    if hasattr(e, "shape"):
        return all(simplify(c) == 0 for c in e)
    return simplify(e) == 0


def check(label, ok):
    if not bool(ok):
        raise AssertionError(f"[FAIL] {label}")
    print(f"  [PASS] {label}")


def jac(v):
    """Lean `grad`: (grad v)_ij = d_j v_i."""
    return Matrix(3, 3, lambda i, j: diff(v[i], coords[j]))


def sym(A):
    return half * (A + A.T)


def skw(A):
    return half * (A - A.T)


def axial(A):
    """Lean `axial`; axial(grad v) = curl v."""
    return Matrix([A[2, 1] - A[1, 2], A[0, 2] - A[2, 0], A[1, 0] - A[0, 1]])


def curl(v):
    return axial(jac(v))


def div(v):
    return jac(v).trace()


def nsq(a):
    return expand((a.T * a)[0])


def generic_symmetric(p):
    s = symbols(f"{p}11 {p}12 {p}13 {p}22 {p}23 {p}33", real=True)
    return Matrix([[s[0], s[1], s[2]], [s[1], s[3], s[4]], [s[2], s[4], s[5]]])


def generic_skew(p):
    a, b, c = symbols(f"{p}1 {p}2 {p}3", real=True)
    return Matrix([[0, -c, b], [c, 0, -a], [-b, a, 0]])


def generic_matrix(p):
    return Matrix(3, 3, symbols(f"{p}11 {p}12 {p}13 {p}21 {p}22 {p}23 {p}31 {p}32 {p}33",
                                real=True))


# Localization: indicator of the critical region (Lean `CriticalRegion`)
Lam, lam = symbols("Lambda lambda", real=True)


def chi_critical(omega):
    return Piecewise((1, Lam**2 < nsq(omega)), (0, True))


# --------------------------------------------------------------------------- A
def helmholtz_uniqueness():
    header("A. Unique sym/skew split of any 3x3 matrix (always true)")
    A = generic_matrix("a")
    S, W = sym(A), skw(A)
    check("A = S(A) + Omega(A)", is_zero(A - S - W))
    check("S(A)^T = S(A)", is_zero(S - S.T))
    check("Omega(A)^T = -Omega(A)", is_zero(W + W.T))
    Sg, Wg = generic_symmetric("s"), generic_skew("w")
    check("uniqueness: sym(S' + W') = S'  (S' symmetric, W' skew)", is_zero(sym(Sg + Wg) - Sg))
    check("uniqueness: skw(S' + W') = W'", is_zero(skw(Sg + Wg) - Wg))
    print("  LOCK A: the Helmholtz split of grad v is unique.      [Lean: helmholtz_unique]")


# --------------------------------------------------------------------------- B
def Z_does_not_see_v():
    header("B. Z == 0 does not mention v")
    v = Matrix([sin(z), 0, 0])
    Sv = sym(jac(v))
    print("  counterexample: X = O = G0 = 0 (so Z == 0) and v = (sin z, 0, 0)")
    check("div v = 0", is_zero(div(v)))
    check("S(v)_13 = cos(z)/2, which is not identically 0 = S(X)",
          is_zero(Sv[0, 2] - cos(z) / 2) and not is_zero(Sv[0, 2]))
    print("  LOCK B: H1 is not a corollary of Z == 0 + div v = 0.")


# --------------------------------------------------------------------------- C
def reconstruction_theorem():
    header("C. Reconstruction postulate (R) => H1 and grad v = G0 (everywhere)")
    G0 = generic_matrix("g")
    # Z == 0 : S(X) = G0 - Omega(O).  By uniqueness its only solution is
    SX, WO = sym(G0), skw(G0)
    check("Z == 0 is solved by S(X) = sym G0, Omega(O) = skw G0", is_zero(SX - (G0 - WO)))
    gradv = SX + WO                                             # (R)
    check("H1: S(v) = S(X)", is_zero(sym(gradv) - SX))
    check("H1: Omega(v) = Omega(O)", is_zero(skw(gradv) - WO))
    check("grad v = G0", is_zero(gradv - G0))
    print("  LOCK C: H1 = uniqueness + (R).   [Lean: H1_of_reconstruction, velocity_gradient_is_G0]")


# --------------------------------------------------------------------------- D
def localized_core():
    header("D. v3: zero-point core localized to CriticalRegion(Lambda)")
    w1, w2, w3 = symbols("omega1 omega2 omega3", real=True)
    omega = Matrix([w1, w2, w3])
    chi = chi_critical(omega)
    print("  chi = 1 on {Lambda^2 < |omega|^2}, 0 elsewhere;  core: chi * (sym G0 - lambda g) = 0")
    check("chi = 1 at a supercritical point (Lambda=1, omega=(2,0,0))",
          chi.subs({Lam: 1, w1: 2, w2: 0, w3: 0}) == 1)
    check("chi = 0 at a subcritical point (Lambda=1, omega=(1/2,0,0))",
          chi.subs({Lam: 1, w1: half, w2: 0, w3: 0}) == 0)

    # D1 inside the core: grad v = G0 with sym G0 = lambda I
    print("\n  D1. inside the core (chi = 1)")
    A = lam * eye(3) + generic_skew("q")                       # grad v = G0, sym G0 = lambda I
    sol = solve(A.trace(), lam)                                 # div v = tr grad v = 3 lambda = 0
    check("div v = 3 lambda = 0  =>  lambda = 0", sol == [0])
    A0 = A.subs(lam, 0)
    check("S(v) = sym(grad v) = 0 on the core", is_zero(sym(A0)))
    check("vortex stretching (grad v) omega = 0 on the core", is_zero(A0 * axial(A0)))
    B = generic_matrix("b")
    check("identity  skw(A) . axial(A) = 0  for every A", is_zero(skw(B) * axial(B)))
    print("  [Lean: lambda_from_incompressibility, strain_vanishes_on_core,")
    print("         skw_mulVec_axial, stretching_vanishes_on_core]")

    # D2 outside the core: no constraint, classical strained flow allowed
    print("\n  D2. outside the core (chi = 0): no constraint")
    A = Matrix([[1, 2, 0], [0, -1, 0], [0, 3, 0]])             # traceless, strained, stretching
    om = axial(A)
    check("example grad v is traceless (div v = 0)", A.trace() == 0)
    check("its vorticity (3, 0, -2) has |omega|^2 = 13 < 16 = Lambda^2 (Lambda = 4)",
          om == Matrix([3, 0, -2]) and nsq(om) == 13)
    c = chi.subs({Lam: 4, w1: om[0], w2: om[1], w3: om[2]})
    check("chi = 0 there, so the core residual chi*(sym G0 - lambda g) vanishes for any lambda",
          c == 0 and is_zero(c * (sym(A) - lam * eye(3))))
    check("strain S(v) != 0 and stretching (grad v) omega = (3,0,0) != 0 are allowed",
          not is_zero(sym(A)) and A * om == Matrix([3, 0, 0]))

    # D3 an explicit flow with vorticity satisfying the localized postulates
    print("\n  D3. Kolmogorov flow v = (sin 2 pi y, 0, 0)   (Lean: x 1 = y)")
    v = Matrix([sin(2 * pi * y), 0, 0])
    J = jac(v)
    om = curl(v)
    check("div v = 0", is_zero(div(v)))
    check("omega = (0, 0, -2 pi cos 2 pi y), not identically zero",
          is_zero(om - Matrix([0, 0, -2 * pi * cos(2 * pi * y)])) and not is_zero(om[2]))
    check("(2 pi)^2 - |omega|^2 = 4 pi^2 sin^2(2 pi y) >= 0, so CriticalRegion(2 pi) is empty",
          is_zero((2 * pi)**2 - nsq(om) - 4 * pi**2 * sin(2 * pi * y)**2))
    print("  velocityState: X = O = v, G0 = grad v, lambda = 0, K = CriticalRegion(2 pi)")
    check("(R): grad v = S(X) + Omega(O)", is_zero(J - sym(J) - skw(J)))
    check("Z == 0: S(X) = G0 - Omega(O)", is_zero(sym(J) - (J - skw(J))))
    print("  [Lean: velocityState_postulates, localized_postulates_admit_vorticity]")

    # D4 the v2 global core as the special case K = T^3
    print("\n  D4. global core (chi == 1 on all of T^3, v2): S(v) == 0 everywhere")
    vf = Matrix([Function(f"v{i}")(x, y, z) for i in (1, 2, 3)])
    e = jac(vf) + jac(vf).T                                     # e_ij = d_j v_i + d_i v_j
    ok = all(expand(2 * diff(vf[k], coords[i], coords[j])
                    - (diff(e[k, j], coords[i]) + diff(e[k, i], coords[j])
                       - diff(e[i, j], coords[k]))) == 0
             for i in range(3) for j in range(3) for k in range(3))
    check("2 d_i d_j v_k = d_i e_kj + d_j e_ki - d_k e_ij  (all 27 index triples)", ok)
    Bm = generic_matrix("m")
    sol = solve(list(Bm * eye(3)), list(Bm))                    # periodic affine: B e_j = 0
    check("Z^3-periodic affine field a + B x forces B = 0", all(val == 0 for val in sol.values()))
    print("  => e == 0 kills all second derivatives: v is affine, hence constant on T^3.")
    print("  [Lean: flat_T3_killing_rigidity, global_core_forces_rigid_motion]")
    print("  LOCK D: the core forces S(v) = 0 only on CriticalRegion; outside it any")
    print("          incompressible flow is admitted.")


# --------------------------------------------------------------------------- E
def rigidity_warning():
    header("E. Do not add a Helmholtz split on the fields X, O themselves")
    print("  Rejected: v = X + O with Omega(X) = 0, S(O) = 0  =>  X = grad phi, O Killing;")
    print("  div O = 0 makes phi harmonic, hence constant on T^3, and the flow a rigid")
    print("  rotation.  Use only (R), a statement about grad v.")
    print("  LOCK E: reconstruction is on grad v, not on the fields X, O.")


def main():
    print("H.U.G.G.E.R CAS -- Proof 1: Helmholtz Bridge (H1)  [v3, localized core]")
    helmholtz_uniqueness()
    Z_does_not_see_v()
    reconstruction_theorem()
    localized_core()
    rigidity_warning()
    header("STATUS")
    print("  Pointwise S(v)=S(X) without (R)        : FALSE   (LOCK B)")
    print("  H1 under reconstruction (R)            : THEOREM [Lean H1_of_reconstruction]")
    print("  Zero-point core imposed on             : CriticalRegion(Lambda) = {Lambda^2 < |omega|^2}")
    print("  On the core: lambda=0, S(v)=0, stretch=0: THEOREM [Lean strain_vanishes_on_core,")
    print("                                                     stretching_vanishes_on_core]")
    print("  Outside the core                       : no constraint (LOCK D2)")
    print("  Flow with omega != 0 admitted          : Kolmogorov, Lambda = 2 pi")
    print("                                           [Lean localized_postulates_admit_vorticity]")
    print("  Global core (v2, K = T^3)              : rigid motion only [Lean global_core_forces_rigid_motion]")
    print("  All checks above passed (assert-guarded).")


if __name__ == "__main__":
    main()
