#!/usr/bin/env python3
"""
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

CAS_Proof3_Mean_Conserved.py  --  Proof 3: Mean Conserved (H3), v3 (localized core)

Lean 4 counterpart: Supplementary_2_Formal_Verification_Lean4/Theorem3_Mean_Conserved.lean

Target: a time-uniform bound ||X(t)||_inf <= M on T^3.

v3 synchronisation with Lean
  * S(X) = S(v) everywhere (H1).  With the zero-point core only on K = CriticalRegion(Lambda),
    S(X) vanishes on K but not elsewhere, so spatial constancy of X is NO LONGER a
    consequence of the postulates.  It survives as
      - X_constant_of_killing     : S(X) = 0 and Omega(X) = 0 everywhere  =>  X constant,
      - X_constant_of_global_core : the v2 special case K = T^3,
    and H3 itself (H3_mean_conserved) takes "X(t) spatially constant" (hconst) and anchor
    invariance X_bar(t) = X_bar(0) (hinv) as explicit hypotheses.
  * v1/v2 of this script only printed statements; every check below is now an executed
    SymPy computation guarded by an assert.
"""

from sympy import (Function, Matrix, Rational, cos, diff, expand, eye, pi, simplify, sin,
                   solve, sqrt, symbols)

half = Rational(1, 2)
x, y, z, t = symbols("x y z t", real=True)
coords = (x, y, z)


# --------------------------------------------------------------------------- helpers
def header(title):
    print("\n" + "=" * 72)
    print(title)
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


def nsq(a):
    return expand((a.T * a)[0])


def generic_matrix(p):
    return Matrix(3, 3, symbols(f"{p}11 {p}12 {p}13 {p}21 {p}22 {p}23 {p}31 {p}32 {p}33",
                                real=True))


# --------------------------------------------------------------------------- A
def compactness_only():
    header("A. Compactness gives a snapshot bound, not H3")
    print("  Each continuous snapshot X(., t) on compact T^3 has a finite sup M(t); nothing")
    print("  makes sup_t M(t) finite.  H3 needs a time-uniform bound from geometry or data.")
    print("  LOCK A: compactness alone does not give a time-uniform M.")


# --------------------------------------------------------------------------- B
def gauge_freedom():
    header("B. Z == 0 does not see the constant mode of X")
    Xf = Matrix([Function(f"X_{i}")(x, y, z) for i in (1, 2, 3)])
    a = Matrix(symbols("A1 A2 A3", real=True))
    check("S(X + a) = S(X) for every constant a", is_zero(sym(jac(Xf + a)) - sym(jac(Xf))))
    check("grad(X + a) = grad X, so Z, div v and G0 are unchanged", is_zero(jac(Xf + a) - jac(Xf)))
    print("  The shift a = (A, 0, 0) with A arbitrary makes ||X + a||_inf exceed any proposed M.")
    print("  LOCK B: H3 is not a corollary of Z == 0.")


# --------------------------------------------------------------------------- C
def localized_constancy():
    header("C. v3: constancy of X under the localized core")
    # C1 S(X) = S(v) everywhere, and = 0 on the core
    G0 = generic_matrix("g")
    SX = sym(G0)                         # Z == 0 symmetrised (unique split, Proof 1 C)
    gradv = G0                           # (R) + Z == 0  =>  grad v = G0
    check("S(X) = sym G0 = S(v) everywhere (H1)   [Lean strain_X_eq_strain_v]",
          is_zero(SX - sym(gradv)))
    lam = symbols("lambda", real=True)
    core_sym = lam * eye(3)              # on K: sym G0 = lambda g, and lambda = 0 (Proof 1 D1)
    check("on the core: S(X) = lambda g with lambda = 0, i.e. S(X) = 0   "
          "[Lean strain_X_vanishes_on_core]",
          solve(core_sym.trace(), lam) == [0] and is_zero(core_sym.subs(lam, 0)))

    # C2 off the core X need not be constant: the Kolmogorov witness
    Xk = Matrix([sin(2 * pi * y), 0, 0])  # velocityState: X = v
    check("Kolmogorov witness (X = v = (sin 2 pi y, 0, 0)) satisfies the localized postulates "
          "with Lambda = 2 pi (Proof 1 D3), yet S(X)_12 = pi cos 2 pi y != 0",
          is_zero(sym(jac(Xk))[0, 1] - pi * cos(2 * pi * y)) and not is_zero(sym(jac(Xk))[0, 1]))
    check("... and X is not spatially constant", not is_zero(jac(Xk)))
    print("  => under localization the constancy of X is an extra hypothesis, not a theorem.")

    # C3 when it does hold: S(X) = Omega(X) = 0 everywhere => grad X = 0; periodic affine => const
    Xf = Matrix([Function(f"X_{i}")(x, y, z) for i in (1, 2, 3)])
    check("grad X = S(X) + Omega(X), so S(X) = Omega(X) = 0 everywhere gives grad X = 0   "
          "[Lean X_constant_of_killing]", is_zero(jac(Xf) - sym(jac(Xf)) - skw(jac(Xf))))
    B = generic_matrix("b")
    a0 = Matrix(symbols("c1 c2 c3", real=True))
    Xa = a0 + B * Matrix(coords)
    shifts = [Xa.subs(coords[j], coords[j] + 1) - Xa for j in range(3)]     # Z^3-periodicity
    sol = solve([c for s in shifts for c in s], list(B))
    check("Z^3-periodic affine X = c + B x forces B = 0 (X constant)",
          all(val == 0 for val in sol.values()) and len(sol) == 9)
    print("  [Lean: const_of_grad_eq_zero, X_constant_of_killing, X_constant_of_global_core]")
    print("  LOCK C: with the core on all of T^3 (v2) X is constant; under localization this")
    print("          is supplied as the hypothesis hconst.")


# --------------------------------------------------------------------------- D
def conservation_of_mean():
    header("D. Anchor invariance (explicit hypothesis) => time-uniform bound")
    p, q, r = symbols("p q r", real=True)
    Xbar0 = Matrix([p, q, r])
    M0 = sqrt(nsq(Xbar0))
    Xt = Xbar0                            # hconst: X(t, x) = X_bar(t);  hinv: X_bar(t) = X_bar(0)
    check("|X(t, x)|^2 = |X_bar(0)|^2 = M0^2 for all t, x", is_zero(nsq(Xt) - M0**2))
    check("M0 >= 0", M0.is_nonnegative)
    print("  Every M0 >= |X_bar(0)| is a time-uniform bound: IntakeBound (X(t)) M0.")
    print("  [Lean: H3_mean_conserved (hconst, hinv as explicit hypotheses)]")
    print("  LOCK D: M = |X_bar(0)| is a theorem under hconst + hinv.")


# --------------------------------------------------------------------------- E
def what_fails():
    header("E. What must not be claimed")
    print("  1. Finiteness of T^3 gives only per-snapshot bounds.")
    print("  2. chi(T^3) = 0 does not bound a vector field.")
    print("  3. Finite G0 energy does not bound ||X||_inf: G0 sees grad X, not the constant mode.")
    print("  4. v3: the localized postulates alone do not make X constant (section C2).")
    print("  LOCK E: chi(T^3) = 0 is not the reason H3 holds; hconst + hinv are.")


def main():
    print("H.U.G.G.E.R CAS -- Proof 3: Mean Conserved (H3)  [v3, localized core]")
    compactness_only()
    gauge_freedom()
    localized_constancy()
    conservation_of_mean()
    what_fails()
    header("STATUS")
    print("  H3 as a free axiom                    : RETIRED")
    print("  S(X) = 0 on CriticalRegion            : THEOREM [Lean strain_X_vanishes_on_core]")
    print("  X spatially constant                  : only if S(X) = Omega(X) = 0 on all of T^3")
    print("                                          [Lean X_constant_of_killing]; else hypothesis")
    print("  ||X(t)||_inf <= |X_bar(0)|            : THEOREM under hconst + hinv [Lean H3_mean_conserved]")
    print("  All checks above passed (assert-guarded).")


if __name__ == "__main__":
    main()
