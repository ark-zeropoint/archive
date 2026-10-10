/-
Copyright (c) 2026 Jung Soo Kim (Ark Project).
SPDX-License-Identifier: CC-BY-NC-4.0
This code is part of the H.U.G.G.E.R + TZT Grand Unified Framework.
See the README.md file in the root directory for full license details.

  File: Theorem0_Basic_Topology.lean
  H.U.G.G.E.R / TZT  —  Layer 0: carrier types, kinematics, TZT state, constants

  v3 (localized zero-point core)
  * The zero-point core `sym 𝒢⁰ = λ g` is imposed only on a core region `K : Set Pt`:
    `ZeroPointCore s lam K` and `TZTPostulates s lam K`.  v2 is the case `K = Set.univ`.
  * `CriticalRegion s Λ = {x | |ω(x)| > Λ}` — the singularity-critical region.
  * `SmoothPeriodic v` = ℤ³-periodic and C² (`ContDiff ℝ 2 v`).
  * §0.6 collects elementary facts about `grad` used throughout.

  Carried over from v2
  * Genuine Mathlib objects: points ℝ³ (universal cover of T³ = ℝ³/ℤ³, fields
    ℤ³-periodic); `(∇v)ᵢⱼ = ∂ⱼvᵢ` from `fderiv`; `dotProduct (⬝ᵥ)`, `crossProduct (⨯₃)`.
  * `λ` is a reserved word in Lean 4, hence `lam`.
  * No global `axiom` anywhere in the project; this file has no `sorry`.
-/

import Mathlib.Analysis.Calculus.ContDiff.Defs
import Mathlib.LinearAlgebra.CrossProduct
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

namespace HUGGER

open scoped Matrix

/-! ## 0.1 Carrier types -/

/-- Euclidean fibre ℝ³ (components indexed by `Fin 3`). -/
abbrev Vec3 : Type := Fin 3 → ℝ

/-- Real 3×3 matrices: velocity-gradient tensors. -/
abbrev Mat3 : Type := Matrix (Fin 3) (Fin 3) ℝ

/-- Points of ℝ³, the universal cover of the flat torus T³ = ℝ³/ℤ³.
A field on T³ is a ℤ³-periodic field on `Pt` (see `Periodic3`). -/
abbrev Pt : Type := Fin 3 → ℝ

/-- ℤ³-periodicity: `f` descends to the flat torus T³ = ℝ³/ℤ³. -/
def Periodic3 {β : Type*} (f : Pt → β) : Prop :=
  ∀ (x : Pt) (n : Fin 3 → ℤ), f (x + fun i => (n i : ℝ)) = f x

/-- Regularity class: ℤ³-periodic and C² (a C² velocity field on T³). -/
structure SmoothPeriodic (v : Pt → Vec3) : Prop where
  periodic : Periodic3 v
  smooth   : ContDiff ℝ 2 v

/-! ## 0.2 Kinematics (Mathlib Fréchet derivative) -/

/-- Velocity gradient `(∇v)ᵢⱼ = ∂ⱼ vᵢ`, read off Mathlib's `fderiv`. -/
noncomputable def grad (v : Pt → Vec3) (x : Pt) : Mat3 :=
  Matrix.of fun i j => fderiv ℝ v x (Pi.single j 1) i

/-- Symmetric part `S(A) = (A + Aᵀ)/2`. -/
noncomputable def sym (A : Mat3) : Mat3 := (1 / 2 : ℝ) • (A + Aᵀ)

/-- Skew part `Ω(A) = (A − Aᵀ)/2`. -/
noncomputable def skw (A : Mat3) : Mat3 := (1 / 2 : ℝ) • (A - Aᵀ)

/-- Strain-rate field `S(v) = sym ∇v`. -/
noncomputable def Strain (v : Pt → Vec3) (x : Pt) : Mat3 := sym (grad v x)

/-- Spin field `Ω(v) = skw ∇v`. -/
noncomputable def Skew (v : Pt → Vec3) (x : Pt) : Mat3 := skw (grad v x)

/-- Divergence `∇·v = tr ∇v`. -/
noncomputable def divergence (v : Pt → Vec3) (x : Pt) : ℝ := Matrix.trace (grad v x)

/-- Axial vector of a matrix; for `A = ∇v` it is the vorticity `∇×v`. -/
def axial (A : Mat3) : Vec3 := ![A 2 1 - A 1 2, A 0 2 - A 2 0, A 1 0 - A 0 1]

/-- Curl `∇×v`. -/
noncomputable def curl (v : Pt → Vec3) (x : Pt) : Vec3 := axial (grad v x)

/-- Squared Euclidean norm `|a|² = a ⬝ᵥ a`. -/
def nsq (a : Vec3) : ℝ := a ⬝ᵥ a

/-! ## 0.3 TZT state and (localized) postulates -/

/-- A TZT snapshot.  The vorticity is not an independent field:
it is *defined* as `ω := ∇×v` (see `TZTState.ω`). -/
structure TZTState where
  /-- intake axis 𝒳 -/
  X  : Pt → Vec3
  /-- twist / rotation field 𝒪 -/
  O  : Pt → Vec3
  /-- velocity -/
  v  : Pt → Vec3
  /-- zero-point metric operator 𝒢⁰ = Λ̂⁰(I) 𝒢 -/
  G0 : Pt → Mat3

/-- Vorticity of a state, `ω = ∇×v`. -/
noncomputable def TZTState.ω (s : TZTState) : Pt → Vec3 := curl s.v

/-- TZT zero-tensor condition `𝒵 ≡ 0`, i.e. `S(𝒳) = 𝒢⁰ − Ω(𝒪)` everywhere. -/
def Z_eq_zero (s : TZTState) : Prop :=
  ∀ x, Strain s.X x = s.G0 x - Skew s.O x

/-- Reconstruction postulate (R): `∇v = S(𝒳) + Ω(𝒪)`.  Constitutive, not derived. -/
def Reconstruction (s : TZTState) : Prop :=
  ∀ x, grad s.v x = Strain s.X x + Skew s.O x

/-- **Localized** zero-point core: on the core region `K`, the symmetric part of `𝒢⁰` is
pure trace, `sym 𝒢⁰(x) = λ(x) g`.  Outside `K` no condition is imposed.
(v2 imposed this on all of T³; that is the special case `K = Set.univ`.) -/
def ZeroPointCore (s : TZTState) (lam : Pt → ℝ) (K : Set Pt) : Prop :=
  ∀ x ∈ K, sym (s.G0 x) = lam x • (1 : Mat3)

/-- Incompressibility `∇·v = 0`. -/
def Incompressible (s : TZTState) : Prop :=
  ∀ x, divergence s.v x = 0

/-- The TZT closure postulates with core region `K`, bundled as *hypotheses* (not global
axioms).  (R), `𝒵 ≡ 0` and `∇·v = 0` hold everywhere; the zero-point core only on `K`. -/
structure TZTPostulates (s : TZTState) (lam : Pt → ℝ) (K : Set Pt) : Prop where
  /-- (R) reconstruction -/
  recon   : Reconstruction s
  /-- 𝒵 ≡ 0 -/
  ztensor : Z_eq_zero s
  /-- zero-point core on `K`: `sym 𝒢⁰ = λ g` -/
  core    : ZeroPointCore s lam K
  /-- `∇·v = 0` -/
  incomp  : Incompressible s

/-- The singularity-critical region at threshold `Λ`: where the vorticity exceeds `Λ`,
`{x | Λ² < |ω(x)|²}`. -/
def CriticalRegion (s : TZTState) (Λ : ℝ) : Set Pt :=
  {x | Λ ^ 2 < nsq (s.ω x)}

/-! ## 0.4 Remainder and intake bound -/

/-- Remainder after the exact index split: `R = (1/4) 𝒳 × (∇×ω)`. -/
noncomputable def Remainder (X ω : Pt → Vec3) (x : Pt) : Vec3 :=
  (1 / 4 : ℝ) • (X x ⨯₃ curl ω x)

/-- Intake-axis structural bound (H3): `0 ≤ M` and `‖𝒳‖_{L∞} ≤ M`
(stated pointwise with the squared Euclidean norm). -/
structure IntakeBound (X : Pt → Vec3) (M : ℝ) : Prop where
  nonneg : 0 ≤ M
  bound  : ∀ x, nsq (X x) ≤ M ^ 2

/-! ## 0.5 Locked constants -/

/-- Viscous absorption constant `K_ν = M² / (16ν)`. -/
noncomputable def Kν (M ν : ℝ) : ℝ := M ^ 2 / (16 * ν)

/-- Young coefficient locked by CAS: `ε* = M² / (32ν)`. -/
noncomputable def epsStar (M ν : ℝ) : ℝ := M ^ 2 / (32 * ν)

/-! ## 0.6 Elementary facts about `grad` -/

/-- Every `u : ℝ³` is `u₀ e₀ + u₁ e₁ + u₂ e₂`. -/
theorem decomp3 (u : Pt) :
    u = u 0 • Pi.single 0 1 + u 1 • Pi.single 1 1 + u 2 • Pi.single 2 1 := by
  funext i
  fin_cases i <;> simp

/-- A continuous linear map on ℝ³ that kills the standard basis is zero. -/
theorem clm_eq_zero_of_single {F : Type*} [NormedAddCommGroup F] [NormedSpace ℝ F]
    (L : Pt →L[ℝ] F) (h : ∀ j : Fin 3, L (Pi.single j 1) = 0) : L = 0 := by
  refine ContinuousLinearMap.ext fun u => ?_
  rw [decomp3 u]
  simp [h]

theorem grad_apply (v : Pt → Vec3) (x : Pt) (i j : Fin 3) :
    grad v x i j = fderiv ℝ v x (Pi.single j 1) i := rfl

/-- `∇v(x) = 0` as a matrix iff the Fréchet derivative vanishes at `x`. -/
theorem grad_eq_zero_iff (v : Pt → Vec3) (x : Pt) : grad v x = 0 ↔ fderiv ℝ v x = 0 := by
  constructor
  · intro h
    refine clm_eq_zero_of_single _ fun j => funext fun i => ?_
    simpa [grad_apply] using congrFun (congrFun h i) j
  · intro h
    ext i j
    simp [grad_apply, h]

theorem axial_zero : axial (0 : Mat3) = 0 := by
  ext i
  fin_cases i <;> simp [axial]

end HUGGER
