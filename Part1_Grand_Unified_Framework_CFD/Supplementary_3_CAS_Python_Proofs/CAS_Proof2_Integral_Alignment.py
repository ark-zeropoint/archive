#!/usr/bin/env python3
"""
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

CAS_Proof2_Integral_Alignment.py  --  Proof 2: Integral Alignment (H2), v3 (localized core)

Lean 4 counterpart: Supplementary_2_Formal_Verification_Lean4/Theorem2_Integral_Alignment.lean

Question: does Z == 0 force  X _| S(omega) = 0, so that the remainder R lives only in the
curl channel?

v3 synchronisation with Lean
  * Lean conventions: (grad w)_ij = d_j w_i, curl = axial(grad).  The remainder of the
    stretching split is  R_i = -1/2 X_j d_i omega_j  (CAS_Proof4, section 2.2).
  * SIGN CORRECTION.  R = -1/2 S(omega) X - 1/4 X x (curl omega).  v1/v2 printed
    "+1/4 X x (curl omega)"; their own code computed -1/4.  Section A now asserts the
    correct sign.  (Lean v3 only uses R = 0 or |R|, so no Lean result depends on it.)
  * The zero-point core holds only on K = CriticalRegion(Lambda).  Hence
      - on int K: grad v is locally constant, curl omega = 0 and R = 0 pointwise
        (Lean remainder_vanishes_on_core_interior, H2_on_core_interior);
      - globally: int <omega, R> = 1/2 int_{T^3 \\ K} omega . S(v) omega, which is no
        longer forced to vanish (v2's "= 0" needed the core on all of T^3).
  * Every check below is an executed SymPy computation guarded by an assert.
"""

from sympy import Function, Matrix, Rational, diff, expand, simplify, symbols

half, quarter = Rational(1, 2), Rational(1, 4)
x, y, z = symbols("x y z", real=True)
coords = (x, y, z)


# --------------------------------------------------------------------------- helpers
def header(t):
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)


def is_zero(e):
    if hasattr(e, "shape"):
        return all(simplify(expand(c)) == 0 for c in e)
    return simplify(expand(e)) == 0


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
    return Matrix([A[2, 1] - A[1, 2], A[0, 2] - A[2, 0], A[1, 0] - A[0, 1]])


def curl(v):
    return axial(jac(v))


def div(v):
    return jac(v).trace()


def vdot(a, b):
    return (a.T * b)[0]


def nsq(a):
    return expand(vdot(a, a))


def remainder(X, w):
    """R_i = -1/2 sum_j X_j d_i w_j  =  -1/2 (grad w)^T X."""
    return -half * jac(w).T * X


X = Matrix([Function(f"X_{i}")(x, y, z) for i in (1, 2, 3)])
w = Matrix([Function(f"omega_{i}")(x, y, z) for i in (1, 2, 3)])


# --------------------------------------------------------------------------- A
def section_identity_split():
    header("A. Exact split of R (no hypothesis)")
    Jw = jac(w)
    R = remainder(X, w)
    Rs = -half * sym(Jw) * X                    # -1/2 X _| S(omega)
    Rk = half * skw(Jw) * X                     # skew part:  -1/2 (Omega^T) X = +1/2 Omega X
    c = curl(w)
    check("R = R_sym + R_skew", is_zero(R - Rs - Rk))
    check("R_sym = -1/2 S(omega) X", is_zero(Rs + half * sym(Jw) * X))
    check("R_skew = -1/4 X x (curl omega)   (= +1/4 (curl omega) x X)",
          is_zero(Rk + quarter * X.cross(c)))
    check("v1/v2 sign '+1/4 X x curl omega' is NOT an identity (difference = -1/2 X x curl omega)",
          not is_zero(Rk - quarter * X.cross(c)) and is_zero((Rk - quarter * X.cross(c))
                                                             + half * X.cross(c)))
    print("  LOCK A: R = -1/2 X _| S(omega) - 1/4 X x (curl omega);  H2 = vanishing of the first term.")


# --------------------------------------------------------------------------- B
def section_Z_does_not_kill_H2():
    header("B. Counterexample: Z == 0 does not force H2")
    M = symbols("M", real=True, nonzero=True)
    f = Function("f")
    Xc = Matrix([M, 0, 0])                      # constant: S(X) = 0, Killing
    wc = Matrix([0, 0, f(x)])
    check("S(X) = 0 for constant X (so Z == 0 holds with O = G0 = 0)", is_zero(sym(jac(Xc))))
    check("div omega = 0", is_zero(div(wc)))
    pair = sym(jac(wc)) * Xc
    check("(X _| S(omega))_3 = (M/2) f'(x), not identically 0",
          is_zero(pair - Matrix([0, 0, M * diff(f(x), x) / 2])))
    print("  LOCK B: H2 is not a corollary of Z == 0 alone.")


# --------------------------------------------------------------------------- C
def section_integral_theorem():
    header("C. Integral identity and its v3 localization")
    R = remainder(X, w)
    SX = sym(jac(X))
    Xw = vdot(X, w)
    lhs = vdot(w, R)
    rhs = (-half * div(Xw * w) + half * Xw * div(w) + half * vdot(w, SX * w))
    check("<omega,R> = -1/2 div((X.omega) omega) + 1/2 (X.omega) div(omega) + 1/2 omega.S(X)omega",
          is_zero(lhs - rhs))
    print("  On T^3 the divergence integrates to 0 and div omega = 0, hence")
    print("      int <omega, R> = 1/2 int omega . S(X) omega,   with S(X) = S(v)  (H1).")

    # v3 localization: chi = indicator of the core region K (1 on K, 0 elsewhere)
    chi = symbols("chi", real=True)
    Sv = symbols("S11 S12 S13 S22 S23 S33", real=True)
    S = Matrix([[Sv[0], Sv[1], Sv[2]], [Sv[1], Sv[3], Sv[4]], [Sv[2], Sv[4], Sv[5]]])
    om = Matrix(symbols("o1 o2 o3", real=True))
    dens = half * vdot(om, (1 - chi) * S * om)  # S(v) = 0 where chi = 1 (Lean strain_vanishes_on_core)
    check("core points (chi = 1): integrand 1/2 omega.S(v)omega = 0",
          is_zero(dens.subs(chi, 1)))
    check("classical points (chi = 0): integrand = 1/2 omega.S(v)omega, unconstrained",
          is_zero(dens.subs(chi, 0) - half * vdot(om, S * om)) and not is_zero(dens.subs(chi, 0)))
    print("  => int <omega,R> = 1/2 int_{T^3 \\ K} omega . S(v) omega : only the classical region")
    print("     contributes, and nothing forces it to vanish.  v2's 'int <omega,R> = 0' required")
    print("     the core on all of T^3 (K = T^3), which admits only rigid motions.")
    print("  LOCK C: integral H2 holds on the core; globally it is NOT forced under localization.")


# --------------------------------------------------------------------------- D
def section_pointwise_under_locks():
    header("D. Pointwise H2 under the three constitutive locks")
    M = symbols("M", real=True, positive=True)
    W = Function("w")
    Xl = Matrix([M, 0, 0])                                       # D1: 0-point frame
    wl = Matrix([W(y, z), 0, 0])                                 # D2: axial vorticity
    pair = sym(jac(wl)) * Xl
    check("D1+D2: X _| S(omega) = (0, (M/2) d_y w, (M/2) d_z w)",
          is_zero(pair - Matrix([0, M * diff(W(y, z), y) / 2, M * diff(W(y, z), z) / 2])))
    wc = symbols("w_c", real=True)
    check("D1+D2+D3 (w constant): X _| S(omega) = 0",
          is_zero(sym(jac(Matrix([wc, 0, 0]))) * Xl))
    print("  LOCK D: pointwise H2 <=> D1+D2+D3 (too rigid: omega constant).")


# --------------------------------------------------------------------------- E
def section_localized_pointwise():
    header("E. v3: pointwise H2 on the interior of the core")
    print("  On an open core set S(v) = 0, so grad v is locally constant (Lean grad_eventually_const):")
    print("  locally v = a + B x with B skew (a rigid motion).")
    a = Matrix(symbols("a1 a2 a3", real=True))
    b1, b2, b3 = symbols("b1 b2 b3", real=True)
    B = Matrix([[0, -b3, b2], [b3, 0, -b1], [-b2, b1, 0]])
    v = a + B * Matrix(coords)
    om = curl(v)
    check("S(v) = 0 for the rigid motion", is_zero(sym(jac(v))))
    check("omega = 2 (b1, b2, b3) is constant", is_zero(om - 2 * Matrix([b1, b2, b3])))
    check("grad omega = 0, hence curl omega = 0 and S(omega) = 0", is_zero(jac(om)))
    check("R = 0 for every intake axis X  =>  <omega, R> = 0 on int K", is_zero(remainder(X, om)))
    print("  [Lean: remainder_vanishes_on_core_interior, H2_on_core_interior]")
    print("  LOCK E: H2 holds pointwise on int CriticalRegion; outside, no alignment is claimed.")


def main():
    print("H.U.G.G.E.R CAS -- Proof 2: Integral Alignment (H2)  [v3, localized core]")
    section_identity_split()
    section_Z_does_not_kill_H2()
    section_integral_theorem()
    section_pointwise_under_locks()
    section_localized_pointwise()
    header("STATUS")
    print("  Remainder split                 : R = -1/2 X_|S(omega) - 1/4 X x curl omega (sign corrected)")
    print("  Pointwise H2 (vector)           : NOT forced by Z == 0 (LOCK B)")
    print("  H2 on int CriticalRegion        : THEOREM, R = 0 [Lean H2_on_core_interior]")
    print("  Integral H2 on all of T^3       : NOT forced under localization; v2 form = global core")
    print("  Lean action                     : H2_integral (global) -> H2_on_core_interior (local)")
    print("  All checks above passed (assert-guarded).")


if __name__ == "__main__":
    main()
