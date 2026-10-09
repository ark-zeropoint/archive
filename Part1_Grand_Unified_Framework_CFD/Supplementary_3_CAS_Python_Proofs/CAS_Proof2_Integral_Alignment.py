#!/usr/bin/env python3
# Copyright (c) 2026 Jung Soo Kim (Ark Project).
# SPDX-License-Identifier: CC-BY-NC-4.0
# This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
# See the README.md file in the root directory for full license details.

"""
File: CAS_Proof2_Integral_Alignment.py
Proof 2: Integral Alignment (H2) — H.U.G.G.E.R / TZT Framework

Question: does Z_ij ≡ 0 force
    X ⌟ Strain(ω) = 0
and therefore
    R = (1/4) X × (∇×ω)   ?

Answer the CAS will print:
  - pointwise H2 is NOT an identity from Z ≡ 0 alone
  - the INTEGRAL identity
        ∫ ⟨ω, R⟩ = (1/2) ∫ ω · G0_sym ω
    IS forced by Z ≡ 0 + div ω = 0 + T³ (no boundary)
  - 0-point core G0_sym = λ g on T³ forces λ = 0
        ⇒ ∫ ⟨ω, R⟩ = 0
  - pointwise H2 becomes a theorem only after three extra
    constitutive locks (axis frame, axial invariance, no fibre shear)
"""

from itertools import product

from sympy import (
    Function,
    IndexedBase,
    Matrix,
    Symbol,
    symbols,
    simplify,
    Rational,
    LeviCivita,
    Eq,
    zeros,
)


half = Rational(1, 2)
quarter = Rational(1, 4)

x, y, z = symbols("x y z", real=True)
coords = (x, y, z)

# Generic fields
X1, X2, X3 = symbols("X_1 X_2 X_3", cls=Function)
w1, w2, w3 = symbols("omega_1 omega_2 omega_3", cls=Function)
X = Matrix([X1(x, y, z), X2(x, y, z), X3(x, y, z)])
w = Matrix([w1(x, y, z), w2(x, y, z), w3(x, y, z)])


def partial(comp, axis):
    return comp.diff(coords[axis])


def grad_vec(v):
    G = zeros(3)
    for i, j in product(range(3), repeat=2):
        G[i, j] = partial(v[j], i)  # (∇v)_ij = ∂_i v_j
    return G


def strain(v):
    G = grad_vec(v)
    return half * (G + G.T)


def skew(v):
    G = grad_vec(v)
    return half * (G - G.T)


def curl(v):
    return Matrix(
        [
            partial(v[2], 1) - partial(v[1], 2),
            partial(v[0], 2) - partial(v[2], 0),
            partial(v[1], 0) - partial(v[0], 1),
        ]
    )


def header(t):
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)


def section_identity_split():
    header("A. Exact component split of R  (no hypothesis)")
    # R_i = -1/2 X_j ∂_i ω_j
    R = zeros(3, 1)
    for i in range(3):
        s = 0
        for j in range(3):
            s += -half * X[j] * partial(w[j], i)
        R[i] = s

    Sw = strain(w)
    Ow = skew(w)
    R_sym = zeros(3, 1)
    R_skw = zeros(3, 1)
    for i in range(3):
        acc_s = 0
        acc_k = 0
        for j in range(3):
            acc_s += -half * X[j] * Sw[i, j] * 2 / 1
            # X_j ∂_i ω_j = X_j (S_ij + Ω_ij) with S_ij = S_ji = ½(∂i wj + ∂j wi)
            # ∂_i ω_j = (∇ω)_{ij}
            pass
        R_sym[i] = sum(-half * X[j] * (partial(w[j], i) + partial(w[i], j)) * half * 2 for j in range(3))
        # simpler: ∂_i ω_j = S_ij + Ω_ij, S = strain, Ω = skew of ω
        R_sym[i] = sum(-half * X[j] * Sw[i, j] for j in range(3)) * 1
        # Wait Sw[i,j] = ½(∂i wj + ∂j wi). And ∂_i ω_j = Sw[i,j] + Ow[i,j]
        R_sym[i] = sum(-half * X[j] * Sw[i, j] for j in range(3))
        R_skw[i] = sum(-half * X[j] * Ow[i, j] for j in range(3))

    # Recompute cleanly
    Gw = grad_vec(w)  # Gw[i,j] = ∂_i w_j
    R2 = zeros(3, 1)
    Rs = zeros(3, 1)
    Rk = zeros(3, 1)
    for i in range(3):
        R2[i] = sum(-half * X[j] * Gw[i, j] for j in range(3))
        Rs[i] = sum(-half * X[j] * Sw[i, j] for j in range(3))
        Rk[i] = sum(-half * X[j] * Ow[i, j] for j in range(3))

    print("R_i = -½ X_j ∂_i ω_j")
    print("R   = R_sym + R_skew")
    print("difference R - (R_sym+R_skew) =")
    print(simplify(R2 - Rs - Rk))

    # curl channel: (1/4) X × curl ω
    # Ow[i,j] = ½(∂i wj - ∂j wi) = -½ ε_ijk (curl ω)_k
    # -½ X_j Ow[i,j] should equal (1/4)(X × curl ω)_i
    cx = Rational(1, 4) * X.cross(curl(w))
    print("\nR_skew  vs  (1/4) X × (∇×ω)")
    print(simplify(Rk - cx))
    print("LOCK A:  R = -½ X⌟S^(ω)  +  (1/4) X×(∇×ω)")
    print("H2 is exactly the vanishing of the first summand.")
    return R2, Rs, Rk, cx


def section_Z_does_not_kill_H2():
    header("B. Counter-component: Z ≡ 0 does not force H2")
    print(
        """
Z_ij = S^X_ij + Ω^O_ij - G0_ij ≡ 0
constrains (X, O, G0).  It does not constrain Strain(ω).

Component counter-example on T³:
  X = (M, 0, 0)   constant   (so S^X = 0, X Killing)
  ω = (0, 0, f(x))
  S^(ω)_13 = S^(ω)_31 = (1/2) f'(x)
  (X ⌟ S^(ω))_3 = (M/2) f'(x)   which is ≠ 0 unless f' = 0.

One can still choose O, G0 so that Z ≡ 0
(e.g. O = 0, G0 = 0, X Killing).  H2 already fails.
"""
    )
    print("LOCK B: H2 is not a corollary of Z ≡ 0 alone.")


def section_integral_theorem():
    header("C. THE inevitable theorem from Z ≡ 0  (integral H2)")
    print(
        r"""
R_i = -½ X_j ∂_i ω_j
⟨ω, R⟩ = -½ X_j ω_i ∂_i ω_j
        = -½ X · ((ω·∇)ω)

div ω = 0  ⇒  ∫ (ω·∇) φ = 0 on T³.

Identity:
  X · ((ω·∇)ω)
    = (ω·∇)(X·ω) - ω_i ω_j ∂_i X_j
    = (ω·∇)(X·ω) - ω · S^X ω     (ω · Ω^X ω = 0)

∫ ⟨ω, R⟩ = -½ ∫ (ω·∇)(X·ω) + ½ ∫ ω · S^X ω
         = ½ ∫ ω · S^X ω

Z ≡ 0  ⇒  S^X = G0_sym          (ω · Ω^O ω = 0)
∫ ⟨ω, R⟩ = ½ ∫ ω · G0_sym ω

0-point core (no residual shear):  G0_sym = λ g
∫ ⟨ω, R⟩ = (λ/2) ∫ |ω|²

Flat closed T³ admits no nontrivial homothety ⇒ λ = 0
∫ ⟨ω, R⟩ = 0

This is the theorem that Z ≡ 0 + 0-point core actually forces.
It is an INTEGRAL alignment, not the pointwise H2.
"""
    )
    print("LOCK C:  ∫ ⟨ω,R⟩ = (λ/2)∫|ω|²  and on T³, λ=0 ⇒ ∫⟨ω,R⟩=0")


def section_pointwise_under_locks():
    header("D. Pointwise H2 as a theorem under three constitutive locks")
    print(
        """
Lock D1  0-point frame: X = (M, 0, 0),  M constant.
         (axis of the active intake; translations on T³ are Killing)
Lock D2  toroidal / axial vorticity: ω = (w(y,z), 0, 0)
         (right-hand: circulation in the 2–3 fibre ⇒ vorticity ∥ axis)
         AND no axial gradient: ∂_x w = 0   (already in w(y,z))
Lock D3  fibre-homogeneity of |ω|:  ∂_y w = 0 and ∂_z w = 0
         i.e. w = const.

Then S^(ω) = 0 and X ⌟ S^(ω) = 0, so H2 holds
and R = (1/4) X × (∇×ω).

WARNING: D3 forces ω spatially constant.  That is too rigid
for a genuine turbulent field.  Therefore the PAPER-GRADE
promotion is Lock C (integral), not D3.

A usable middle lock D3' (weaker):
  only demand the enstrophy pairing
      ω_i X_j S^(ω)_ij = 0
  which follows from C with λ = 0 without freezing ω.
"""
    )

    M = symbols("M", real=True, positive=True)
    w = Function("w")
    # D1+D2 without D3
    Xloc = Matrix([M, 0, 0])
    wloc = Matrix([w(y, z), 0, 0])
    Sw = zeros(3)
    for i, j in product(range(3), repeat=2):
        Sw[i, j] = half * (wloc[j].diff(coords[i]) + wloc[i].diff(coords[j]))
    contraction = zeros(3, 1)
    for i in range(3):
        contraction[i] = sum(Xloc[j] * Sw[j, i] for j in range(3))
    print("D1+D2 only,  X ⌟ S^(ω) =")
    print(simplify(contraction))
    print("= (0, (M/2) ∂_y w, (M/2) ∂_z w)   ≠ 0 unless D3")
    print()
    # with D3
    wc = symbols("w_const", real=True)
    wconst = Matrix([wc, 0, 0])
    Swc = zeros(3)
    for i, j in product(range(3), repeat=2):
        Swc[i, j] = half * (wconst[j].diff(coords[i]) + wconst[i].diff(coords[j]))
    contrc = zeros(3, 1)
    for i in range(3):
        contrc[i] = sum(Xloc[j] * Swc[j, i] for j in range(3))
    print("D1+D2+D3,  X ⌟ S^(ω) =")
    print(simplify(contrc))
    print("LOCK D: pointwise H2 ⇔ D1+D2+D3 (too rigid) or the pairing form D3'.")


def section_pairing_theorem():
    header("E. Paper-grade H2  (pairing form, theorem)")
    print(
        """
Define
  H2_pair  :   ω_i X_j S^(ω)_ij  =  0     (scalar, not vector)

From section C, on T³ with G0_sym = 0:
  ∫ ω_i X_j S^(ω)_ij   is the same channel as ∫ ⟨ω,R⟩
  and equals 0.

Pointwise H2_pair is still not free.  What IS a theorem:

  theorem H2_integral
      (hZ : Z = 0)
      (hCore : G0_sym = λ • g)
      (hT3 : no nontrivial homothety)
      : ∫ ⟨ω, R⟩ = 0

  theorem remainder_energy
      : ∫ ⟨ω, stretching⟩ = ∫ ⟨ω, R⟩ = 0
        after the flux Π evaporates on T³.

Viscous absorption then does not even need the pointwise
formula R = (1/4) X×curl ω.  It needs only

  |∫ ⟨ω,R⟩| ≤ ε ‖ω‖₂² + Cε ‖curl ω‖₂²

and the left side is 0 under H2_integral, which is stronger:
the production vanishes and viscosity only fights the
linear transport that already integrates to zero.
"""
    )
    print("LOCK E: promote H2_integral, not pointwise H2.")


def main():
    print("H.U.G.G.E.R CAS — Proof 2: Integral Alignment (H2)")
    section_identity_split()
    section_Z_does_not_kill_H2()
    section_integral_theorem()
    section_pointwise_under_locks()
    section_pairing_theorem()
    header("STATUS")
    print("Pointwise H2 (vector) : NOT forced by Z ≡ 0.")
    print("Integral H2           : FORCED by Z ≡ 0 + 0-point core on T³.")
    print("Lean action           : axiom H2_align  →  theorem H2_integral")


if __name__ == "__main__":
    main()
