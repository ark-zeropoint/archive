/-
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

  File: Theorem2_Integral_Alignment.lean
  H2 promotion + rigidity analysis.

  CAS: Supplementary_3_CAS_Python_Proofs/CAS_Proof2_Integral_Alignment.py

  Pointwise vector H2  X ⌟ Strain(ω) = 0  is NOT a corollary of Z ≡ 0.

  v3
  * §2.3 local rigidity: if S(v) = 0 near x (v of class C²) then ∇²v(x) = 0.
    Proof: ∂ₖ∂ⱼvᵢ is symmetric in (k, j) (Schwarz, Mathlib `isSymmSndFDerivAt`) and
    antisymmetric in (j, i) (differentiate S(v) = 0), hence zero.
  * §2.4 global rigidity on flat T³ (v2's `sorry`, now proved): a C² ℤ³-periodic
    field with S(v) ≡ 0 is a constant translation (∇v ≡ 0).
  * §2.5 H2 localized: on the interior of the core region the remainder channel is
    shut, R = 0, hence ⟨ω, R⟩ = 0 there.  Outside the core nothing is claimed.
  * §2.6 the v2 global statement is the special case K = Set.univ.
  No axiom, no sorry.
-/

import Theorem1_Helmholtz_Bridge
import Mathlib.Analysis.Calculus.ContDiff.Comp
import Mathlib.Analysis.Calculus.FDeriv.Symmetric
import Mathlib.Analysis.Calculus.MeanValue

namespace HUGGER

open scoped Matrix Topology
open Filter

/-! ## 2.1 `λ = 0` on a nonempty core (v1's constant-`λ` form) -/

theorem T3_no_homothety (s : TZTState) (lam : ℝ) (K : Set Pt)
    (hP : TZTPostulates s (fun _ => lam) K) (hK : K.Nonempty) : lam = 0 := by
  obtain ⟨x, hx⟩ := hK
  exact lambda_from_incompressibility s _ K hP hx

/-! ## 2.2 Consequences of C² regularity -/

theorem differentiable_of_C2 {v : Pt → Vec3} (hv : ContDiff ℝ 2 v) : Differentiable ℝ v :=
  hv.differentiable (by norm_num)

theorem contDiff_fderiv_of_C2 {v : Pt → Vec3} (hv : ContDiff ℝ 2 v) :
    ContDiff ℝ 1 (fderiv ℝ v) :=
  hv.fderiv_right (m := 1) (by norm_num)

theorem differentiable_fderiv_of_C2 {v : Pt → Vec3} (hv : ContDiff ℝ 2 v) :
    Differentiable ℝ (fderiv ℝ v) :=
  (contDiff_fderiv_of_C2 hv).differentiable (by norm_num)

/-! ## 2.3 Local rigidity: vanishing strain kills the Hessian -/

/-- Evaluation functional `L ↦ (L eⱼ)ᵢ`. -/
noncomputable def evalE (i j : Fin 3) : (Pt →L[ℝ] Vec3) →L[ℝ] ℝ :=
  (ContinuousLinearMap.proj i).comp (ContinuousLinearMap.apply ℝ Vec3 (Pi.single j 1))

@[simp] theorem evalE_apply (i j : Fin 3) (L : Pt →L[ℝ] Vec3) :
    evalE i j L = L (Pi.single j 1) i := rfl

/-- If `v` is C² and its strain rate vanishes on a neighbourhood of `x`,
then the second derivative of `v` vanishes at `x`. -/
theorem hessian_eq_zero_of_eventually_strain_zero (v : Pt → Vec3) (hv : ContDiff ℝ 2 v)
    (x : Pt) (hS : ∀ᶠ y in 𝓝 x, Strain v y = 0) : fderiv ℝ (fderiv ℝ v) x = 0 := by
  have hdiff : DifferentiableAt ℝ (fderiv ℝ v) x := differentiable_fderiv_of_C2 hv x
  have hsymm : IsSymmSndFDerivAt ℝ v x := hv.contDiffAt.isSymmSndFDerivAt (by simp)
  set D2 := fderiv ℝ (fderiv ℝ v) x with hD2
  -- antisymmetry in the last two slots: differentiate `∂ⱼvᵢ + ∂ᵢvⱼ = 0`
  have hanti : ∀ k j i : Fin 3,
      D2 (Pi.single k 1) (Pi.single j 1) i + D2 (Pi.single k 1) (Pi.single i 1) j = 0 := by
    intro k j i
    have h1 : HasFDerivAt (fun y => evalE i j (fderiv ℝ v y) + evalE j i (fderiv ℝ v y))
        ((evalE i j).comp D2 + (evalE j i).comp D2) x :=
      ((evalE i j).hasFDerivAt.comp x hdiff.hasFDerivAt).add
        ((evalE j i).hasFDerivAt.comp x hdiff.hasFDerivAt)
    have h0 : (fun y => evalE i j (fderiv ℝ v y) + evalE j i (fderiv ℝ v y)) =ᶠ[𝓝 x]
        fun _ => (0 : ℝ) := by
      filter_upwards [hS] with y hy
      have h := congrFun (congrFun hy i) j
      simp only [Strain, sym, Matrix.smul_apply, Matrix.add_apply, Matrix.transpose_apply,
        Matrix.zero_apply, smul_eq_mul, grad_apply] at h
      simp only [evalE_apply]
      linarith
    have h2 : HasFDerivAt (fun y => evalE i j (fderiv ℝ v y) + evalE j i (fderiv ℝ v y))
        (0 : Pt →L[ℝ] ℝ) x :=
      (hasFDerivAt_const (0 : ℝ) x).congr_of_eventuallyEq h0
    have h3 := congrArg (fun L : Pt →L[ℝ] ℝ => L (Pi.single k 1)) (h1.unique h2)
    simpa using h3
  -- symmetry in the first two slots (Schwarz)
  have hsym : ∀ k j i : Fin 3,
      D2 (Pi.single k 1) (Pi.single j 1) i = D2 (Pi.single j 1) (Pi.single k 1) i :=
    fun k j i => congrFun (hsymm (Pi.single k 1) (Pi.single j 1)) i
  -- a 3-tensor symmetric in (k, j) and antisymmetric in (j, i) vanishes
  have hT : ∀ k j i : Fin 3, D2 (Pi.single k 1) (Pi.single j 1) i = 0 := by
    intro k j i
    have e1 := hanti k j i
    have e2 := hsym k i j
    have e3 := hanti i k j
    have e4 := hsym i j k
    have e5 := hanti j i k
    have e6 := hsym j k i
    linarith
  refine clm_eq_zero_of_single D2 fun k => clm_eq_zero_of_single _ fun j => ?_
  funext i
  exact hT k j i

/-- On an open set where the strain rate vanishes, `∇v` is locally constant. -/
theorem grad_eventually_const (v : Pt → Vec3) (hv : ContDiff ℝ 2 v) {U : Set Pt}
    (hU : IsOpen U) (hS : ∀ y ∈ U, Strain v y = 0) {x : Pt} (hx : x ∈ U) :
    ∀ᶠ y in 𝓝 x, grad v y = grad v x := by
  obtain ⟨r, hr, hball⟩ := Metric.isOpen_iff.mp hU x hx
  have hD : Differentiable ℝ (fderiv ℝ v) := differentiable_fderiv_of_C2 hv
  have hH : Set.EqOn (fderiv ℝ (fderiv ℝ v)) 0 (Metric.ball x r) := fun y hy =>
    hessian_eq_zero_of_eventually_strain_zero v hv y
      (eventually_of_mem (hU.mem_nhds (hball hy)) hS)
  filter_upwards [Metric.ball_mem_nhds x hr] with y hy
  have h := Metric.isOpen_ball.is_const_of_fderiv_eq_zero (convex_ball x r).isPreconnected
    hD.differentiableOn hH hy (Metric.mem_ball_self hr)
  ext i j
  simp only [grad_apply, h]

/-! ## 2.4 Global rigidity of the flat torus (proved; was `sorry` in v2) -/

/-- **Rigidity of flat T³.**  A C² ℤ³-periodic field with vanishing strain rate
(a Killing field of flat T³) has vanishing gradient: it is a constant translation. -/
theorem flat_T3_killing_rigidity (v : Pt → Vec3) (hv : SmoothPeriodic v)
    (hS : ∀ x, Strain v x = 0) (x : Pt) : grad v x = 0 := by
  have hD : Differentiable ℝ (fderiv ℝ v) := differentiable_fderiv_of_C2 hv.smooth
  have hH : ∀ y, fderiv ℝ (fderiv ℝ v) y = 0 := fun y =>
    hessian_eq_zero_of_eventually_strain_zero v hv.smooth y (Eventually.of_forall hS)
  -- the derivative is a constant linear map `L`
  have hconst : ∀ y, fderiv ℝ v y = fderiv ℝ v 0 := fun y =>
    is_const_of_fderiv_eq_zero hD hH y 0
  set L := fderiv ℝ v 0 with hL
  -- `v` is affine: `v y = v 0 + L y`
  have haff : ∀ y, v y = v 0 + L y := by
    have hg : Differentiable ℝ (fun y => v y - L y) :=
      (differentiable_of_C2 hv.smooth).sub L.differentiable
    have hg' : ∀ y, fderiv ℝ (fun y => v y - L y) y = 0 := by
      intro y
      rw [fderiv_fun_sub (differentiable_of_C2 hv.smooth y) L.differentiableAt, L.fderiv,
        hconst y, sub_self]
    intro y
    have h := is_const_of_fderiv_eq_zero hg hg' y 0
    rw [map_zero, sub_zero] at h
    exact sub_eq_iff_eq_add.mp h
  -- ℤ³-periodicity forces `L = 0`
  have hLe : ∀ j : Fin 3, L (Pi.single j 1) = 0 := by
    intro j
    have hcast : ((0 : Pt) + fun i => ((Pi.single j (1 : ℤ) : Fin 3 → ℤ) i : ℝ)) =
        Pi.single j 1 := by
      funext i
      by_cases h : i = j
      · subst h
        simp
      · simp [h]
    have hp : v (Pi.single j 1) = v 0 := by
      have := hv.periodic 0 (Pi.single j 1)
      rwa [hcast] at this
    rw [haff (Pi.single j 1)] at hp
    simpa using hp
  rw [grad_eq_zero_iff, hconst x, clm_eq_zero_of_single L hLe]

/-! ## 2.5 H2 on the core region -/

/-- On the interior of the core region the remainder channel is shut: `R = 0`. -/
theorem remainder_vanishes_on_core_interior (s : TZTState) (lam : Pt → ℝ) (K : Set Pt)
    (hP : TZTPostulates s lam K) (hv : SmoothPeriodic s.v) {x : Pt}
    (hx : x ∈ interior K) : Remainder s.X s.ω x = 0 := by
  have hev := grad_eventually_const s.v hv.smooth isOpen_interior
    (fun y hy => strain_vanishes_on_core s lam K hP (interior_subset hy)) hx
  have hω : s.ω =ᶠ[𝓝 x] fun _ => s.ω x := by
    filter_upwards [hev] with y hy
    show axial (grad s.v y) = axial (grad s.v x)
    rw [hy]
  have hfd : fderiv ℝ s.ω x = 0 := by
    rw [hω.fderiv_eq]
    exact fderiv_const_apply _
  have hcurl : curl s.ω x = 0 := by
    show axial (grad s.ω x) = 0
    rw [(grad_eq_zero_iff _ _).mpr hfd, axial_zero]
  simp [Remainder, hcurl]

/-- **H2 (localized)**: the production density `⟨ω, R⟩` vanishes on the interior of the
core region.  Outside the core no alignment is claimed. -/
theorem H2_on_core_interior (s : TZTState) (lam : Pt → ℝ) (K : Set Pt)
    (hP : TZTPostulates s lam K) (hv : SmoothPeriodic s.v) {x : Pt}
    (hx : x ∈ interior K) : s.ω x ⬝ᵥ Remainder s.X s.ω x = 0 := by
  rw [remainder_vanishes_on_core_interior s lam K hP hv hx, dotProduct_zero]

/-! ## 2.6 The v2 global core as a special case (`K = Set.univ`) -/

/-- With the zero-point core imposed everywhere, the vorticity vanishes identically. -/
theorem vorticity_vanishes_of_global_core (s : TZTState) (lam : Pt → ℝ)
    (hP : TZTPostulates s lam Set.univ) (hv : SmoothPeriodic s.v) (x : Pt) : s.ω x = 0 := by
  have hg := flat_T3_killing_rigidity s.v hv
    (fun y => strain_vanishes_on_core s lam Set.univ hP (Set.mem_univ y)) x
  show axial (grad s.v x) = 0
  rw [hg, axial_zero]

end HUGGER
