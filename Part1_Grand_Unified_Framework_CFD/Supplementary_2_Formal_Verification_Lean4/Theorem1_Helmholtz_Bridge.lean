/-
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

  File: Theorem1_Helmholtz_Bridge.lean
  H1 promotion.

  CAS: Supplementary_3_CAS_Python_Proofs/CAS_Proof1_Helmholtz_Bridge.py

  Z ≡ 0 does not mention v, so S(v) = S(X) is not a corollary of Z ≡ 0 alone.
  Reconstruction postulate (R):  ∇v = S(X) + Ω(O).
  Helmholtz (sym/skew) uniqueness ⇒ H1 is a theorem (everywhere).
  Z ≡ 0 + (R) ⇒ ∇v = G0 (everywhere).

  v3 (localized core): on the core region K,
    div v = 0  ⇒  λ = 0  ⇒  S(v) = 0  ⇒  vortex stretching (ω·∇)v = (∇v)ω = 0.
  Outside K nothing is forced, so flows with vorticity are admitted.
  All statements are kernel-checked theorems (no axiom, no sorry).
-/

import Theorem0_Basic_Topology

namespace HUGGER

open scoped Matrix

/-! ## 1.1 Pointwise sym/skew algebra -/

theorem sym_add_skw (A : Mat3) : sym A + skw A = A := by
  ext i j
  simp only [sym, skw, Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply,
    Matrix.transpose_apply, smul_eq_mul]
  ring

theorem sym_transpose (A : Mat3) : (sym A)ᵀ = sym A := by
  ext i j
  simp only [sym, Matrix.transpose_apply, Matrix.smul_apply, Matrix.add_apply, smul_eq_mul]
  ring

theorem skw_transpose (A : Mat3) : (skw A)ᵀ = -skw A := by
  ext i j
  simp only [skw, Matrix.transpose_apply, Matrix.neg_apply, Matrix.smul_apply,
    Matrix.sub_apply, smul_eq_mul]
  ring

theorem sym_sym (A : Mat3) : sym (sym A) = sym A := by
  ext i j
  simp only [sym, Matrix.add_apply, Matrix.smul_apply, Matrix.transpose_apply, smul_eq_mul]
  ring

theorem sym_skw (A : Mat3) : sym (skw A) = 0 := by
  ext i j
  simp only [sym, skw, Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply,
    Matrix.transpose_apply, Matrix.zero_apply, smul_eq_mul]
  ring

theorem sym_sub (A B : Mat3) : sym (A - B) = sym A - sym B := by
  ext i j
  simp only [sym, Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply,
    Matrix.transpose_apply, smul_eq_mul]
  ring

theorem trace_sym (A : Mat3) : Matrix.trace (sym A) = Matrix.trace A := by
  rw [sym, Matrix.trace_smul, Matrix.trace_add, Matrix.trace_transpose, smul_eq_mul]
  ring

/-- Unique Helmholtz split: `A = S + W`, `Sᵀ = S`, `Wᵀ = −W` ⇒ `S = sym A`, `W = skw A`. -/
theorem helmholtz_unique (A S W : Mat3) (h : A = S + W) (hS : Sᵀ = S) (hW : Wᵀ = -W) :
    S = sym A ∧ W = skw A := by
  subst h
  have hS' : ∀ i j, S j i = S i j := fun i j => by
    simpa only [Matrix.transpose_apply] using congrFun (congrFun hS i) j
  have hW' : ∀ i j, W j i = -W i j := fun i j => by
    simpa only [Matrix.transpose_apply, Matrix.neg_apply] using congrFun (congrFun hW i) j
  constructor
  · ext i j
    simp only [sym, Matrix.add_apply, Matrix.smul_apply, Matrix.transpose_apply, smul_eq_mul]
    rw [hS' i j, hW' i j]
    ring
  · ext i j
    simp only [skw, Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply,
      Matrix.transpose_apply, smul_eq_mul]
    rw [hS' i j, hW' i j]
    ring

/-- The spin part annihilates the axial vector: `Ω(A) · axial(A) = 0`
(for `A = ∇v` this is `½ ω × ω = 0`). -/
theorem skw_mulVec_axial (A : Mat3) : skw A *ᵥ axial A = 0 := by
  ext i
  fin_cases i <;>
    simp [skw, axial, Matrix.mulVec, dotProduct, Fin.sum_univ_three] <;> ring

/-! ## 1.2 H1 and ∇v = 𝒢⁰ (everywhere) -/

/-- **H1** (a theorem): under (R) the NS strain/spin coincide with the TZT fields. -/
theorem H1_of_reconstruction (s : TZTState) (hR : Reconstruction s) (x : Pt) :
    Strain s.v x = Strain s.X x ∧ Skew s.v x = Skew s.O x := by
  obtain ⟨h1, h2⟩ := helmholtz_unique (grad s.v x) (Strain s.X x) (Skew s.O x) (hR x)
    (sym_transpose _) (skw_transpose _)
  exact ⟨h1.symm, h2.symm⟩

/-- `𝒵 ≡ 0` + (R) ⇒ `∇v = 𝒢⁰`. -/
theorem velocity_gradient_is_G0 (s : TZTState) (hR : Reconstruction s) (hZ : Z_eq_zero s)
    (x : Pt) : grad s.v x = s.G0 x := by
  rw [hR x, hZ x, sub_add_cancel]

/-! ## 1.3 Consequences on the core region `K` -/

/-- On the core: `div v = tr ∇v = tr 𝒢⁰ = tr (sym 𝒢⁰) = 3λ`, so incompressibility
forces `λ = 0` at every core point. -/
theorem lambda_from_incompressibility (s : TZTState) (lam : Pt → ℝ) (K : Set Pt)
    (hP : TZTPostulates s lam K) {x : Pt} (hx : x ∈ K) : lam x = 0 := by
  have h2 : Matrix.trace (grad s.v x) = 3 * lam x := by
    rw [velocity_gradient_is_G0 s hP.recon hP.ztensor x, ← trace_sym, hP.core x hx,
      Matrix.trace_smul, Matrix.trace_one, Fintype.card_fin, smul_eq_mul]
    push_cast
    ring
  have h3 : 3 * lam x = 0 := by
    rw [← h2]
    exact hP.incomp x
  linarith

/-- On the core the strain rate of `v` vanishes: `S(v)(x) = 0` for `x ∈ K`. -/
theorem strain_vanishes_on_core (s : TZTState) (lam : Pt → ℝ) (K : Set Pt)
    (hP : TZTPostulates s lam K) {x : Pt} (hx : x ∈ K) : Strain s.v x = 0 := by
  show sym (grad s.v x) = 0
  rw [velocity_gradient_is_G0 s hP.recon hP.ztensor x, hP.core x hx,
    lambda_from_incompressibility s lam K hP hx, zero_smul]

/-- **Vortex stretching vanishes in the core**: `(ω·∇)v = (∇v) ω = 0` on `K`.
(The only source term of the enstrophy equation is switched off inside the core.) -/
theorem stretching_vanishes_on_core (s : TZTState) (lam : Pt → ℝ) (K : Set Pt)
    (hP : TZTPostulates s lam K) {x : Pt} (hx : x ∈ K) : grad s.v x *ᵥ s.ω x = 0 := by
  have hS : sym (grad s.v x) = 0 := strain_vanishes_on_core s lam K hP hx
  show grad s.v x *ᵥ axial (grad s.v x) = 0
  calc grad s.v x *ᵥ axial (grad s.v x)
      = (sym (grad s.v x) + skw (grad s.v x)) *ᵥ axial (grad s.v x) := by
        rw [sym_add_skw]
    _ = skw (grad s.v x) *ᵥ axial (grad s.v x) := by rw [hS, zero_add]
    _ = 0 := skw_mulVec_axial _

end HUGGER
