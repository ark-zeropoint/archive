/-
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

  File: Theorem5_Consistency_Audit.lean

  5.1 Consistency.  Every hypothesis bundle used in files 0–4 is satisfiable;
      explicit witnesses are constructed below.
  5.2 (v3) The localized postulates admit flows WITH vorticity:
      * `velocityState_postulates` — every C² periodic divergence-free flow whose
        vorticity stays below Λ satisfies `TZTPostulates` with core region
        `K = CriticalRegion Λ` (taking 𝒳 = 𝒪 = v, 𝒢⁰ = ∇v);
      * `localized_postulates_admit_vorticity` — the Kolmogorov shear flow
        v = (sin 2πx₁, 0, 0), with ω = (0, 0, −2π cos 2πx₁) ≢ 0, is such a model.
  5.3 Audit of the global core (v2 behaviour, K = Set.univ): ∇v ≡ 0 and ω ≡ 0.
  No axiom, no sorry.
-/

import Theorem2_Integral_Alignment
import Theorem3_Mean_Conserved
import Theorem4_Viscous_Absorption
import Mathlib.Analysis.Calculus.ContDiff.Operations
import Mathlib.Analysis.Calculus.FDeriv.Mul
import Mathlib.Analysis.Calculus.FDeriv.Prod
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv

namespace HUGGER

open scoped Matrix

/-! ## 5.1 Consistency witnesses -/

theorem grad_const (c : Vec3) (x : Pt) : grad (fun _ => c) x = 0 := by
  ext i j
  simp [grad]

theorem sym_zero : sym (0 : Mat3) = 0 := by
  simp [sym]

theorem skw_zero : skw (0 : Mat3) = 0 := by
  simp [skw]

/-- Constant-velocity state: `v ≡ c`, `𝒳 ≡ a`, `𝒪 ≡ 0`, `𝒢⁰ ≡ 0`. -/
noncomputable def constState (a c : Vec3) : TZTState where
  X  := fun _ => a
  O  := fun _ => 0
  v  := fun _ => c
  G0 := fun _ => 0

/-- The constant states satisfy every TZT postulate, even with the core everywhere. -/
theorem constState_postulates (a c : Vec3) :
    TZTPostulates (constState a c) (fun _ => 0) Set.univ where
  recon x := by
    simp only [constState, Strain, Skew, grad_const, sym_zero, skw_zero, add_zero]
  ztensor x := by
    simp only [constState, Strain, Skew, grad_const, sym_zero, skw_zero, sub_zero]
  core x _ := by
    simp only [constState, sym_zero, zero_smul]
  incomp x := by
    simp only [constState, divergence, grad_const, Matrix.trace_zero]

theorem constState_smooth (a c : Vec3) : SmoothPeriodic (constState a c).v where
  periodic _ _ := rfl
  smooth := contDiff_const

theorem constState_intake (a c : Vec3) (M : ℝ) (hM : 0 ≤ M) (ha : nsq a ≤ M ^ 2) :
    IntakeBound (constState a c).X M :=
  ⟨hM, fun _ => ha⟩

/-- The Step-3 budget hypotheses are consistent (the zero budget satisfies them). -/
noncomputable def zeroBudget (ν M T : ℝ) : EnstrophyBudget ν M T where
  E := fun _ => 0
  E' := fun _ => 0
  D := fun _ => 0
  P := fun _ => 0
  deriv_E _ _ := hasDerivAt_const _ _
  E_nonneg _ := le_rfl
  D_nonneg _ := le_rfl
  balance _ _ := by simp
  holder _ _ := by simp

/-! ## 5.2 The localized postulates admit flows with vorticity -/

/-- The TZT state carried by a velocity field: `𝒳 = 𝒪 = v`, `𝒢⁰ = ∇v`. -/
noncomputable def velocityState (v : Pt → Vec3) : TZTState where
  X  := v
  O  := v
  v  := v
  G0 := grad v

/-- **Converse.** Every divergence-free flow whose vorticity stays below `Λ` satisfies the
TZT postulates with core region `CriticalRegion Λ` (the core condition is then vacuous). -/
theorem velocityState_postulates (v : Pt → Vec3) (hdiv : ∀ x, divergence v x = 0) (Λ : ℝ)
    (hb : ∀ x, nsq (curl v x) ≤ Λ ^ 2) :
    TZTPostulates (velocityState v) (fun _ => 0) (CriticalRegion (velocityState v) Λ) where
  recon x := (sym_add_skw (grad v x)).symm
  ztensor x := eq_sub_of_add_eq (sym_add_skw (grad v x))
  core x hx := absurd hx (not_lt.mpr (hb x))
  incomp := hdiv

/-- Kolmogorov shear flow `v(x) = (sin 2πx₁, 0, 0)`. -/
noncomputable def kolmogorov : Pt → Vec3 :=
  fun x => Real.sin (2 * Real.pi * x 1) • (Pi.single 0 1 : Vec3)

theorem grad_kolmogorov (x : Pt) (i j : Fin 3) :
    grad kolmogorov x i j = Real.cos (2 * Real.pi * x 1) *
      (2 * Real.pi * (Pi.single j 1 : Pt) 1) * (Pi.single 0 1 : Vec3) i := by
  have h1 := (hasFDerivAt_apply (𝕜 := ℝ) (1 : Fin 3) x).const_mul (2 * Real.pi)
  have h2 : HasFDerivAt kolmogorov _ x := h1.sin.smul_const (Pi.single 0 1 : Vec3)
  rw [grad_apply, h2.fderiv]
  simp only [ContinuousLinearMap.smulRight_apply, smul_apply,
    ContinuousLinearMap.proj_apply, Pi.smul_apply, smul_eq_mul]

theorem kolmogorov_div (x : Pt) : divergence kolmogorov x = 0 := by
  simp [divergence, Matrix.trace, Fin.sum_univ_three, grad_kolmogorov]

theorem kolmogorov_curl (x : Pt) :
    curl kolmogorov x = ![0, 0, -(2 * Real.pi * Real.cos (2 * Real.pi * x 1))] := by
  ext i
  fin_cases i
  · simp [curl, axial, grad_kolmogorov]
  · simp [curl, axial, grad_kolmogorov]
  · simp [curl, axial, grad_kolmogorov]
    ring

theorem kolmogorov_bound (x : Pt) : nsq (curl kolmogorov x) ≤ (2 * Real.pi) ^ 2 := by
  have hc := Real.cos_sq_le_one (2 * Real.pi * x 1)
  have h := mul_le_mul_of_nonneg_left hc (sq_nonneg (2 * Real.pi))
  simp only [kolmogorov_curl, nsq, dotProduct, Fin.sum_univ_three]
  simp
  nlinarith

theorem kolmogorov_curl_ne_zero : curl kolmogorov 0 ≠ 0 := by
  intro h
  have h2 := congrFun h 2
  simp [kolmogorov_curl, Real.pi_ne_zero] at h2

theorem kolmogorov_smoothPeriodic : SmoothPeriodic kolmogorov where
  periodic x n := by
    simp only [kolmogorov, Pi.add_apply]
    rw [show 2 * Real.pi * (x 1 + (n 1 : ℝ)) = 2 * Real.pi * x 1 + (n 1 : ℤ) * (2 * Real.pi) by
      ring, Real.sin_add_int_mul_two_pi]
  smooth :=
    (Real.contDiff_sin.comp (contDiff_const.mul (contDiff_apply ℝ ℝ 1))).smul contDiff_const

/-- **The localized postulates admit flows with vorticity.**  There is a C² periodic state
satisfying `TZTPostulates` (with a core region) whose vorticity is not identically zero. -/
theorem localized_postulates_admit_vorticity :
    ∃ (s : TZTState) (lam : Pt → ℝ) (K : Set Pt),
      TZTPostulates s lam K ∧ SmoothPeriodic s.v ∧ ∃ x, s.ω x ≠ 0 :=
  ⟨velocityState kolmogorov, fun _ => 0,
    CriticalRegion (velocityState kolmogorov) (2 * Real.pi),
    velocityState_postulates kolmogorov kolmogorov_div _ kolmogorov_bound,
    kolmogorov_smoothPeriodic, 0, kolmogorov_curl_ne_zero⟩

/-! ## 5.3 Audit of the global core (v2 behaviour, `K = Set.univ`) -/

/-- With the zero-point core imposed on all of T³, every regular state is a constant
translation: `∇v ≡ 0` and `ω ≡ 0`.  (This is why v3 localizes the core.) -/
theorem global_core_forces_rigid_motion (s : TZTState) (lam : Pt → ℝ)
    (hP : TZTPostulates s lam Set.univ) (hv : SmoothPeriodic s.v) :
    (∀ x, grad s.v x = 0) ∧ (∀ x, s.ω x = 0) :=
  ⟨flat_T3_killing_rigidity s.v hv
      (fun y => strain_vanishes_on_core s lam Set.univ hP (Set.mem_univ y)),
    vorticity_vanishes_of_global_core s lam hP hv⟩

end HUGGER
