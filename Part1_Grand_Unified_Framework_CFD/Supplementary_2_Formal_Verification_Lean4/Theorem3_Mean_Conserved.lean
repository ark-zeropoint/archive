/-
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

  File: Theorem3_Mean_Conserved.lean
  H3 promotion.

  CAS: Supplementary_3_CAS_Python_Proofs/CAS_Proof3_Mean_Conserved.py

  Compactness of T³ bounds each snapshot, not the orbit in time.
  Z ≡ 0 is invariant under X ↦ X + a, so it cannot cap ‖X‖∞ by itself.

  v3 (localized core)
  * S(X) = S(v) everywhere (H1).  On the core region K it vanishes; outside K it
    need not, so X is spatially constant only when its strain and spin vanish on
    all of T³ (`X_constant_of_killing`; v2's global-core case is
    `X_constant_of_global_core`).
  * `const_of_grad_eq_zero` (v2's `sorry`) is now proved from Mathlib's
    `is_const_of_fderiv_eq_zero`.
  * H3 itself is unchanged: spatial constancy + anchor invariance X̄(t) = X̄(0)
    (explicit TZT dynamical hypothesis) ⇒ ‖X(t)‖∞ ≤ M₀ for every M₀ ≥ |X̄(0)|.
  No axiom, no sorry.
-/

import Theorem1_Helmholtz_Bridge
import Mathlib.Analysis.Calculus.MeanValue

namespace HUGGER

open scoped Matrix

/-- `S(𝒳) = S(v)` everywhere (from H1). -/
theorem strain_X_eq_strain_v (s : TZTState) (hR : Reconstruction s) (x : Pt) :
    Strain s.X x = Strain s.v x :=
  (H1_of_reconstruction s hR x).1.symm

/-- On the core, symmetrising `𝒵 ≡ 0` gives `S(𝒳) = λ g`. -/
theorem strain_X_eq_core (s : TZTState) (lam : Pt → ℝ) (K : Set Pt) (hZ : Z_eq_zero s)
    (hC : ZeroPointCore s lam K) {x : Pt} (hx : x ∈ K) :
    Strain s.X x = lam x • (1 : Mat3) := by
  have h := congrArg sym (hZ x)
  simp only [Strain, Skew, sym_sym, sym_sub, sym_skw, sub_zero] at h
  show sym (grad s.X x) = lam x • (1 : Mat3)
  rw [h, hC x hx]

/-- On the core the intake axis is strain-free. -/
theorem strain_X_vanishes_on_core (s : TZTState) (lam : Pt → ℝ) (K : Set Pt)
    (hP : TZTPostulates s lam K) {x : Pt} (hx : x ∈ K) : Strain s.X x = 0 := by
  rw [strain_X_eq_strain_v s hP.recon x]
  exact strain_vanishes_on_core s lam K hP hx

/-- A field with identically vanishing gradient on ℝ³ is constant
(proved from Mathlib's `is_const_of_fderiv_eq_zero`). -/
theorem const_of_grad_eq_zero (X : Pt → Vec3) (hX : Differentiable ℝ X)
    (h : ∀ x, grad X x = 0) (x y : Pt) : X x = X y :=
  is_const_of_fderiv_eq_zero hX (fun z => (grad_eq_zero_iff X z).mp (h z)) x y

/-- Strain-free and spin-free everywhere ⇒ `𝒳` is spatially constant. -/
theorem X_constant_of_killing (X : Pt → Vec3) (hX : Differentiable ℝ X)
    (hS : ∀ x, Strain X x = 0) (hW : ∀ x, Skew X x = 0) (x y : Pt) : X x = X y := by
  refine const_of_grad_eq_zero X hX (fun z => ?_) x y
  rw [← sym_add_skw (grad X z)]
  show Strain X z + Skew X z = 0
  rw [hS z, hW z, add_zero]

/-- v2's statement as the global-core special case (`K = Set.univ`):
`λ = 0` everywhere and `Ω(𝒳) = 0` ⇒ `𝒳` is spatially constant. -/
theorem X_constant_of_global_core (s : TZTState) (lam : Pt → ℝ)
    (hP : TZTPostulates s lam Set.univ) (hKill : ∀ x, Skew s.X x = 0)
    (hX : Differentiable ℝ s.X) (x y : Pt) : s.X x = s.X y :=
  X_constant_of_killing s.X hX
    (fun z => strain_X_vanishes_on_core s lam Set.univ hP (Set.mem_univ z)) hKill x y

/-- **H3** (mean conservation ⇒ uniform intake bound).
Hypotheses: `𝒳(t)` is spatially constant with value `X̄(t)` and the topological anchor
is invariant, `X̄(t) = X̄(0)` (v1 `axiom mean_invariant`, now an explicit TZT dynamical
hypothesis).  Then every `M₀ ≥ |X̄(0)|` bounds `𝒳(t)` for all `t`. -/
theorem H3_mean_conserved (s : ℝ → TZTState) (Xbar : ℝ → Vec3)
    (hconst : ∀ t x, (s t).X x = Xbar t) (hinv : ∀ t, Xbar t = Xbar 0)
    (M0 : ℝ) (hM0 : 0 ≤ M0) (h0 : nsq (Xbar 0) ≤ M0 ^ 2) :
    ∀ t, IntakeBound (s t).X M0 := by
  intro t
  refine ⟨hM0, fun x => ?_⟩
  rw [hconst t x, hinv t]
  exact h0

end HUGGER
