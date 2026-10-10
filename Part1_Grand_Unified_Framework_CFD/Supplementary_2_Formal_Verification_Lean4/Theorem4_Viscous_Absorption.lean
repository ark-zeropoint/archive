/-
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

  File: Theorem4_Viscous_Absorption.lean
  Step 3 — Young absorption, viscous absorption, Grönwall (L² envelope)

  CAS counterpart: Supplementary_3_CAS_Python_Proofs/CAS_Proof4_Viscous_Absorption.py
  Locked constants:
    epsStar M ν = M^2 / (32 * ν)
    Kν      M ν = M^2 / (16 * ν)

  v2 changes
  * `viscosity_swallows_half` now assumes `M ≠ 0` (it is false at M = 0, because
    ε* = 0 there and Lean's convention is x / 0 = 0).
  * v1 axioms H1_bridge / H2_align / H3_bound are gone: H1–H3 are theorems in
    files 1–3.  v1 `H3_bound` (∀ M, IntakeBound X M) was inconsistent (take M = −1).
  * The PDE-level facts consumed by Step 3 (enstrophy identity, Hölder bound) are
    the explicit fields of `EnstrophyBudget`, not global axioms.
  * Young and viscous absorption are fully proved.
  v3: the scalar Grönwall lemma (v2's `sorry`) is proved from Mathlib
      (`antitoneOn_of_deriv_nonpos` applied to t ↦ f(t)·e^{−Kt}).  No sorry.

  Step 3 stops at the L² envelope.  The Hˢ/Sobolev route to the BKM integral is
  not formalized; the BKM bound under the localized core is `BKM_sup_bound`
  (Theorem 6), obtained directly from the vorticity ceiling.
-/

import Theorem0_Basic_Topology
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.SpecialFunctions.ExpDeriv

namespace HUGGER

/-! ## 4.0 Constants (kernel-checked) -/

theorem epsStar_eq (M ν : ℝ) : epsStar M ν = M ^ 2 / (32 * ν) := rfl

theorem Kν_eq (M ν : ℝ) : Kν M ν = M ^ 2 / (16 * ν) := rfl

theorem Kν_from_epsStar (M ν : ℝ) : Kν M ν = 2 * epsStar M ν := by
  unfold Kν epsStar
  ring

theorem epsStar_pos (M ν : ℝ) (hM : M ≠ 0) (hν : 0 < ν) : 0 < epsStar M ν := by
  unfold epsStar
  positivity

/-- `ν − M²/(64 ε*) = ν/2`: viscosity swallows exactly half (needs `M ≠ 0`, `ν ≠ 0`). -/
theorem viscosity_swallows_half (M ν : ℝ) (hM : M ≠ 0) (hν : ν ≠ 0) :
    ν - M ^ 2 / (64 * epsStar M ν) = ν / 2 := by
  unfold epsStar
  field_simp
  ring

/-! ## 4.1 Young absorption (scalar, fully proved) -/

/-- Scalar Young step: if `P² ≤ (M²/16)·E·D` with `E, D ≥ 0`, then for every `ε > 0`,
`|P| ≤ ε E + (M²/(64ε)) D`. -/
theorem young_enstrophy (P E D M ε : ℝ) (hε : 0 < ε) (hE : 0 ≤ E) (hD : 0 ≤ D)
    (hP : P ^ 2 ≤ M ^ 2 / 16 * E * D) :
    |P| ≤ ε * E + M ^ 2 / (64 * ε) * D := by
  have hb : 0 ≤ M ^ 2 / (64 * ε) * D := by positivity
  have ha : 0 ≤ ε * E := by positivity
  have hε0 : ε ≠ 0 := hε.ne'
  have hab : ε * E * (M ^ 2 / (64 * ε) * D) = M ^ 2 / 64 * E * D := by
    field_simp
  have key : P ^ 2 ≤ (ε * E + M ^ 2 / (64 * ε) * D) ^ 2 := by
    nlinarith [sq_nonneg (ε * E - M ^ 2 / (64 * ε) * D)]
  have habs := sq_le_sq.mp key
  rwa [abs_of_nonneg (add_nonneg ha hb)] at habs

/-! ## 4.2 Enstrophy budget: PDE-level facts as explicit hypotheses -/

/-- Scalar enstrophy budget of one solution on `[0, T]`.
These are the PDE facts that Step 3 consumes.  In this skeleton they are
**assumptions** (fields), not consequences of the vorticity equation:
* `deriv_E`  : `E(t) = ‖ω(t)‖₂²` is differentiable with derivative `E'(t)`;
* `balance`  : `½ E' + ν D = P` — enstrophy identity after the flux `∇·Π`
               integrates out on T³ (`D = ‖∇×ω‖₂²`, `P = ∫⟨ω, R⟩`);
* `holder`   : `P² ≤ (M²/16) E D` — Cauchy–Schwarz with `‖𝒳‖∞ ≤ M` and
               `R = ¼ 𝒳 × ∇×ω`. -/
structure EnstrophyBudget (ν M T : ℝ) where
  E  : ℝ → ℝ
  E' : ℝ → ℝ
  D  : ℝ → ℝ
  P  : ℝ → ℝ
  deriv_E  : ∀ t ∈ Set.Icc 0 T, HasDerivAt E (E' t) t
  E_nonneg : ∀ t, 0 ≤ E t
  D_nonneg : ∀ t, 0 ≤ D t
  balance  : ∀ t ∈ Set.Icc 0 T, 1 / 2 * E' t + ν * D t = P t
  holder   : ∀ t ∈ Set.Icc 0 T, P t ^ 2 ≤ M ^ 2 / 16 * E t * D t

/-! ## 4.3 Viscous absorption (fully proved) -/

/-- With `ε = ε*`:  `E' ≤ K_ν E − ν D` on `[0, T]`. -/
theorem viscous_absorption {ν M T : ℝ} (B : EnstrophyBudget ν M T) (hν : 0 < ν)
    (hM : M ≠ 0) : ∀ t ∈ Set.Icc 0 T, B.E' t ≤ Kν M ν * B.E t - ν * B.D t := by
  intro t ht
  have hY := young_enstrophy (B.P t) (B.E t) (B.D t) M (epsStar M ν)
    (epsStar_pos M ν hM hν) (B.E_nonneg t) (B.D_nonneg t) (B.holder t ht)
  have hhalf : M ^ 2 / (64 * epsStar M ν) = ν / 2 := by
    have := viscosity_swallows_half M ν hM hν.ne'
    linarith
  rw [hhalf] at hY
  have hbal := B.balance t ht
  have hle : B.P t ≤ |B.P t| := le_abs_self _
  rw [Kν_from_epsStar]
  linarith

/-! ## 4.4 Grönwall (proved from Mathlib) and the L² envelope -/

/-- Scalar Grönwall inequality: `f' ≤ K f` on `[0, T]` ⇒ `f(t) ≤ f(0) e^{Kt}`.
Proof: `g(t) = f(t) e^{−Kt}` has `g' = (f' − K f) e^{−Kt} ≤ 0`, so `g` is antitone. -/
theorem gronwall_scalar (K T : ℝ) (f f' : ℝ → ℝ)
    (hf : ∀ t ∈ Set.Icc 0 T, HasDerivAt f (f' t) t)
    (hle : ∀ t ∈ Set.Icc 0 T, f' t ≤ K * f t) :
    ∀ t ∈ Set.Icc 0 T, f t ≤ f 0 * Real.exp (K * t) := by
  set g : ℝ → ℝ := fun t => f t * Real.exp (-(K * t)) with hg
  have hgd : ∀ t ∈ Set.Icc 0 T, HasDerivAt g
      (f' t * Real.exp (-(K * t)) + f t * (Real.exp (-(K * t)) * -(K * 1))) t := by
    intro t ht
    exact (hf t ht).mul (((hasDerivAt_id' t).const_mul K).neg.exp)
  have hanti : AntitoneOn g (Set.Icc 0 T) := by
    apply antitoneOn_of_deriv_nonpos (convex_Icc 0 T)
    · exact fun t ht => (hgd t ht).continuousAt.continuousWithinAt
    · intro t ht
      rw [interior_Icc] at ht
      exact (hgd t (Set.Ioo_subset_Icc_self ht)).differentiableAt.differentiableWithinAt
    · intro t ht
      rw [interior_Icc] at ht
      have ht' := Set.Ioo_subset_Icc_self ht
      rw [(hgd t ht').deriv]
      have he := Real.exp_pos (-(K * t))
      have hl := hle t ht'
      nlinarith
  intro t ht
  have h0 : (0 : ℝ) ∈ Set.Icc 0 T := ⟨le_rfl, ht.1.trans ht.2⟩
  have hgt := hanti h0 ht ht.1
  simp only [hg, mul_zero, neg_zero, Real.exp_zero, mul_one] at hgt
  have key : f t = f t * Real.exp (-(K * t)) * Real.exp (K * t) := by
    rw [mul_assoc, ← Real.exp_add, neg_add_cancel, Real.exp_zero, mul_one]
  rw [key]
  exact mul_le_mul_of_nonneg_right hgt (Real.exp_pos (K * t)).le

/-- **Step 3 conclusion** (L² envelope): `‖ω(t)‖₂² ≤ ‖ω(0)‖₂² · exp(K_ν t)` on `[0, T]`. -/
theorem enstrophy_gronwall_L2 {ν M T : ℝ} (B : EnstrophyBudget ν M T) (hν : 0 < ν)
    (hM : M ≠ 0) : ∀ t ∈ Set.Icc 0 T, B.E t ≤ B.E 0 * Real.exp (Kν M ν * t) := by
  refine gronwall_scalar (Kν M ν) T B.E B.E' B.deriv_E (fun t ht => ?_)
  have h := viscous_absorption B hν hM t ht
  have hνD : 0 ≤ ν * B.D t := mul_nonneg hν.le (B.D_nonneg t)
  linarith

end HUGGER
