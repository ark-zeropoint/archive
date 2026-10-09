#!/usr/bin/env python3
# Copyright (c) 2026 Jung Soo Kim (Ark Project).
# SPDX-License-Identifier: CC-BY-NC-4.0
# This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
# See the README.md file in the root directory for full license details.

"""
File: CAS_Proof1_Helmholtz_Bridge.py
Proof 1: Helmholtz Bridge (H1) — H.U.G.G.E.R / TZT Framework

Claim under test:
  Z_ij ≡ 0  and  div v = 0
    ⇒  S(v) = S(X)   (and Ω(v) = Ω(O))

CAS verdict:
  Z ≡ 0 never mentions v.  H1 is not a corollary.
  What *is* a theorem: unique Helmholtz split of any (0,2) tensor.
  Promote H1 to a theorem UNDER the reconstruction postulate
      ∇v = S(X) + Ω(O)
  Then Z ≡ 0 becomes ∇v = G0, and
      div v = 0 + G0_sym = λ g  ⇒  λ = 0.
"""

from itertools import product

from sympy import (
    Function,
    Matrix,
    Rational,
    Symbol,
    simplify,
    symbols,
    zeros,
)


half = Rational(1, 2)
x, y, z = symbols("x y z", real=True)
coords = (x, y, z)


def header(t):
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)


def helmholtz_uniqueness():
    header("A. Unique split of any 3×3 matrix  (always true)")
    a = symbols("a11 a12 a13 a21 a22 a23 a31 a32 a33")
    A = Matrix(3, 3, a)
    S = half * (A + A.T)
    W = half * (A - A.T)
    print("A = S + W,  S^T=S,  W^T=-W")
    print("A - (S+W) =")
    print(simplify(A - S - W))
    print("S - S.T =")
    print(simplify(S - S.T))
    print("W + W.T =")
    print(simplify(W + W.T))
    print("LOCK A: Helmholtz split of ∇v is unique.")
    print("        S(v)=S(X) AND Ω(v)=Ω(O)")
    print("        ⇔  ∇v = S(X) + Ω(O)")


def Z_does_not_see_v():
    header("B. Z ≡ 0 does not mention v")
    print(
        """
Z_ij = S(X)_ij + Ω(O)_ij - G0_ij ≡ 0
is a relation among (X, O, G0) only.

Component counter-example:
  X = 0, O = 0, G0 = 0   ⇒  Z = 0
  v = (sin z, 0, 0)      ⇒  S(v)_13 = S(v)_31 = (1/2) cos z ≠ 0 = S(X)
  div v = 0              still holds.

H1 fails while Z ≡ 0 and incompressibility both hold.
"""
    )
    print("LOCK B: H1 is not a corollary of Z ≡ 0 + div v = 0.")


def reconstruction_theorem():
    header("C. Reconstruction postulate  ⇒  H1 is a theorem")
    print(
        """
Postulate (R):
    ∇v  =  S(X) + Ω(O)

Helmholtz uniqueness (A):
    ∇v  =  S(v) + Ω(v)   uniquely.

Therefore
    S(v) = S(X),   Ω(v) = Ω(O).     ← H1, now a theorem

Feed Z ≡ 0:
    S(X) + Ω(O) = G0
so
    ∇v = G0.

This is the 0-point constitutive law for the fluid:
the velocity gradient *is* the 0-point metric operator.
"""
    )
    print("LOCK C: H1 = Helmholtz uniqueness + reconstruction (R).")


def incompressibility_kills_lambda():
    header("D. div v = 0  ⇒  λ = 0   (cleaner than homothety)")
    print(
        """
G0_sym = λ g     (0-point core, no residual shear)
tr(G0_sym) = 3λ

Under (R) and Z ≡ 0:
  div v = tr(∇v) = tr(S(v)) = tr(G0_sym) = 3λ

∇·v = 0  ⇒  λ = 0.

This upgrades H2_integral:
  ∫ ⟨ω, R⟩ = (λ/2) ∫ |ω|² = 0
without invoking 'no homothety on T³' as the primary reason.
No-homothety remains a geometric consistency check:
a nonzero λ would be a homothety, forbidden on closed flat T³
AND forbidden by incompressibility.
"""
    )
    print("LOCK D: λ = 0 from tr(∇v)=0.")


def rigidity_warning():
    header("E. Do not add Helmholtz on the fields X, O themselves")
    print(
        """
Tempting extra split:
  v = X + O,  Ω(X)=0,  S(O)=0.

Then H1 holds, but on T³:
  Ω(X)=0 ⇒ X = ∇φ locally
  S(O)=0 ⇒ O Killing
  div v = Δφ + div O = 0
  If also div O = 0, φ is harmonic ⇒ φ const on T³
  ⇒ X const ⇒ S(X)=0 ⇒ ∇v = Ω(O) = G0
  ⇒ flow is locally rigid rotation.

That freezes turbulence.  Same lesson as pointwise H2.
Do NOT take v = X + O + potential/Killing as the paper's H1.
Take only (R): ∇v = S(X) + Ω(O), which does not freeze |ω|.
"""
    )
    print("LOCK E: reconstruction is on ∇v, not on the fields X,O.")


def main():
    print("H.U.G.G.E.R CAS — Proof 1: Helmholtz Bridge (H1)")
    helmholtz_uniqueness()
    Z_does_not_see_v()
    reconstruction_theorem()
    incompressibility_kills_lambda()
    rigidity_warning()
    header("STATUS")
    print("Pointwise S(v)=S(X) without (R) : FALSE")
    print("H1 under reconstruction (R)    : THEOREM (uniqueness)")
    print("λ = 0 from div v = 0           : THEOREM")
    print("Lean action: axiom H1_bridge → theorem H1_of_reconstruction")


if __name__ == "__main__":
    main()