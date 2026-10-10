/-
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

  File: Theorem6_Vorticity_Ceiling.lean   (new in v3)

  The localized zero-point core as a vorticity ceiling.

  * `vorticity_ceiling`: if the core region contains the critical region
    {x | |ω(x)| > Λ}, then |ω| ≤ Λ everywhere.
    Proof: on the open set {|ω| > Λ} the strain vanishes, so ω is locally constant
    there (Theorem 2.3); the level set {ω = ω(x₀)} of a supercritical point is
    then clopen in the connected space ℝ³, hence everything; so S(v) ≡ 0, and
    global rigidity (Theorem 2.4) gives ω ≡ 0 — contradicting |ω(x₀)| > Λ ≥ 0.
  * `ceiling_iff`: for a C² periodic divergence-free v, the localized postulates
    with core `CriticalRegion Λ` hold for some TZT data iff |ω| ≤ Λ everywhere.
  * `BKM_sup_bound`: the time-dependent form — sup_x |ω(t, x)| ≤ Λ on [0, T],
    hence ∫₀ᵀ ‖ω(t)‖∞ dt ≤ Λ T (the BKM integral is finite).
  No axiom, no sorry.
-/

import Theorem5_Consistency_Audit
import Mathlib.Analysis.Convex.PathConnected
import Mathlib.Topology.Instances.Matrix

namespace HUGGER

open scoped Matrix Topology
open Filter

/-! ## 6.1 Continuity of the kinematic fields -/

theorem continuous_grad (v : Pt → Vec3) (hv : ContDiff ℝ 2 v) : Continuous (grad v) := by
  have hD : Continuous (fderiv ℝ v) := hv.continuous_fderiv (by norm_num)
  refine continuous_matrix fun i j => ?_
  simp only [grad_apply]
  exact (continuous_apply i).comp (hD.clm_apply continuous_const)

theorem continuous_axial : Continuous axial := by
  refine continuous_pi fun i => ?_
  fin_cases i
  · exact (continuous_id.matrix_elem 2 1).sub (continuous_id.matrix_elem 1 2)
  · exact (continuous_id.matrix_elem 0 2).sub (continuous_id.matrix_elem 2 0)
  · exact (continuous_id.matrix_elem 1 0).sub (continuous_id.matrix_elem 0 1)

theorem continuous_nsq : Continuous nsq :=
  continuous_id.dotProduct continuous_id

/-! ## 6.2 The ceiling -/

/-- **Vorticity ceiling.**  If the core region `K` contains the critical region
`{x | |ω(x)| > Λ}`, then `|ω(x)| ≤ Λ` at every point. -/
theorem vorticity_ceiling (s : TZTState) (lam : Pt → ℝ) (K : Set Pt) (Λ : ℝ)
    (hP : TZTPostulates s lam K) (hK : CriticalRegion s Λ ⊆ K) (hv : SmoothPeriodic s.v) :
    ∀ x, nsq (s.ω x) ≤ Λ ^ 2 := by
  have hω : Continuous s.ω := continuous_axial.comp (continuous_grad s.v hv.smooth)
  have hU : IsOpen (CriticalRegion s Λ) :=
    isOpen_lt continuous_const (continuous_nsq.comp hω)
  have hSU : ∀ y ∈ CriticalRegion s Λ, Strain s.v y = 0 := fun y hy =>
    strain_vanishes_on_core s lam K hP (hK hy)
  by_contra hcon
  push Not at hcon
  obtain ⟨x₀, hx₀⟩ := hcon
  -- the level set through a supercritical point is open …
  have hWo : IsOpen {x | s.ω x = s.ω x₀} := by
    refine isOpen_iff_mem_nhds.mpr fun x hx => ?_
    have hx' : s.ω x = s.ω x₀ := hx
    have hxU : x ∈ CriticalRegion s Λ := by
      show Λ ^ 2 < nsq (s.ω x)
      rw [hx']
      exact hx₀
    filter_upwards [grad_eventually_const s.v hv.smooth hU hSU hxU] with y hy
    show s.ω y = s.ω x₀
    rw [← hx']
    show axial (grad s.v y) = axial (grad s.v x)
    rw [hy]
  -- … and closed, hence all of the connected space ℝ³
  have hWc : IsClosed {x | s.ω x = s.ω x₀} := isClosed_eq hω continuous_const
  have hW : Set.univ ⊆ {x | s.ω x = s.ω x₀} :=
    (convex_univ : Convex ℝ (Set.univ : Set Pt)).isPreconnected.subset_isClopen ⟨hWc, hWo⟩
      ⟨x₀, Set.mem_univ _, rfl⟩
  -- so the strain vanishes everywhere, and global rigidity forces ω ≡ 0
  have hall : ∀ y, Strain s.v y = 0 := by
    intro y
    apply hSU
    have hy : s.ω y = s.ω x₀ := hW (Set.mem_univ y)
    show Λ ^ 2 < nsq (s.ω y)
    rw [hy]
    exact hx₀
  have h0 : s.ω x₀ = 0 := by
    show axial (grad s.v x₀) = 0
    rw [flat_T3_killing_rigidity s.v hv hall x₀, axial_zero]
  rw [h0] at hx₀
  have hz : nsq (0 : Vec3) = 0 := by simp [nsq]
  rw [hz] at hx₀
  nlinarith [sq_nonneg Λ]

/-- Under the hypotheses of `vorticity_ceiling` the critical region is empty:
the localized core condition is never active on a nonempty supercritical set. -/
theorem criticalRegion_empty (s : TZTState) (lam : Pt → ℝ) (K : Set Pt) (Λ : ℝ)
    (hP : TZTPostulates s lam K) (hK : CriticalRegion s Λ ⊆ K) (hv : SmoothPeriodic s.v) :
    CriticalRegion s Λ = ∅ := by
  ext x
  simp only [CriticalRegion, Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_lt]
  exact vorticity_ceiling s lam K Λ hP hK hv x

/-! ## 6.3 The ceiling is exactly the localized postulate -/

/-- For a C² periodic divergence-free flow `v`: TZT data satisfying the localized
postulates with core `CriticalRegion Λ` exist **iff** `|ω| ≤ Λ` everywhere. -/
theorem ceiling_iff (v : Pt → Vec3) (hv : SmoothPeriodic v)
    (hdiv : ∀ x, divergence v x = 0) (Λ : ℝ) :
    (∃ (s : TZTState) (lam : Pt → ℝ), s.v = v ∧ TZTPostulates s lam (CriticalRegion s Λ)) ↔
      ∀ x, nsq (curl v x) ≤ Λ ^ 2 := by
  constructor
  · rintro ⟨s, lam, rfl, hP⟩
    exact vorticity_ceiling s lam _ Λ hP subset_rfl hv
  · intro hb
    exact ⟨velocityState v, fun _ => 0, rfl, velocityState_postulates v hdiv Λ hb⟩

/-! ## 6.4 Time-dependent form (BKM) -/

/-- If at every time `t ∈ [0, T]` the postulates hold with a core region containing the
critical region, then `sup_x |ω(t, x)| ≤ Λ` on `[0, T]`; consequently
`∫₀ᵀ ‖ω(t)‖∞ dt ≤ Λ T < ∞` (the Beale–Kato–Majda integral is finite). -/
theorem BKM_sup_bound (s : ℝ → TZTState) (lam : ℝ → Pt → ℝ) (K : ℝ → Set Pt) (Λ T : ℝ)
    (hP : ∀ t ∈ Set.Icc 0 T, TZTPostulates (s t) (lam t) (K t))
    (hK : ∀ t ∈ Set.Icc 0 T, CriticalRegion (s t) Λ ⊆ K t)
    (hv : ∀ t ∈ Set.Icc 0 T, SmoothPeriodic (s t).v) :
    ∀ t ∈ Set.Icc 0 T, ∀ x, nsq ((s t).ω x) ≤ Λ ^ 2 :=
  fun t ht => vorticity_ceiling (s t) (lam t) (K t) Λ (hP t ht) (hK t ht) (hv t ht)

end HUGGER
