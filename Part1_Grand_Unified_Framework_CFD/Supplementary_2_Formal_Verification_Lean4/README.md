# H.U.G.G.E.R + TZT — Lean 4 formal skeleton (Supplementary 2, v3)

Lean **v4.34.1** + Mathlib **v4.34.1**. `lake build` finishes with **0 errors and 0 warnings**.
The project contains **no `axiom` and no `sorry`**. `lake env lean AxiomAudit.lean` reports
that each of the 35 audited results depends only on Lean's standard axioms
`propext`, `Classical.choice` and `Quot.sound`.

```bash
lake exe cache get              # download prebuilt Mathlib
lake build                      # Theorem0 … Theorem6
lake env lean AxiomAudit.lean   # axiom dependencies of every main result
```

## v3: the zero-point core is localized

In v2 the core condition was imposed on all of T³, which forced every admissible flow to be a
constant translation. In v3 it holds only on a core region `K`:

* `ZeroPointCore s lam K`: `sym 𝒢⁰(x) = λ(x) g` for `x ∈ K` only.
* `TZTPostulates s lam K`: (R), `𝒵 ≡ 0` and `∇·v = 0` everywhere; the zero-point core on `K`.
  v2 is the special case `K = Set.univ`.
* `CriticalRegion s Λ = {x | |ω(x)| > Λ}`: the singularity-critical region.

| result | statement | file |
|---|---|---|
| core algebra | on `K`: `λ = 0`, `S(v) = 0`, and vortex stretching `(ω·∇)v = (∇v)ω = 0` | 1 |
| local rigidity | if `S(v) = 0` on a neighbourhood of `x`, then `∇²v(x) = 0`; so on an open set with `S(v) = 0`, `∇v` is locally constant | 2 |
| H2 (localized) | on `interior K`: `R = 0`, hence `⟨ω, R⟩ = 0` | 2 |
| flows with vorticity are admitted | every C² periodic divergence-free `v` with `|ω| ≤ Λ` satisfies the postulates with core `CriticalRegion Λ` (`velocityState_postulates`) | 5 |
| explicit example | the Kolmogorov flow `v = (sin 2πx₁, 0, 0)` with `ω = (0, 0, −2π cos 2πx₁) ≢ 0` (`localized_postulates_admit_vorticity`) | 5 |
| vorticity ceiling | if `CriticalRegion Λ ⊆ K`, then `|ω| ≤ Λ` everywhere (`vorticity_ceiling`) and the critical region is empty (`criticalRegion_empty`) | 6 |
| equivalence | for C² periodic divergence-free `v`: TZT data with core `CriticalRegion Λ` exist **iff** `|ω| ≤ Λ` (`ceiling_iff`) | 6 |
| BKM form | if the postulates hold on `[0, T]` with `CriticalRegion Λ ⊆ K(t)`, then `sup_x |ω(t, x)| ≤ Λ`, so `∫₀ᵀ ‖ω‖∞ dt ≤ ΛT` (`BKM_sup_bound`) | 6 |
| global core (v2) | `K = Set.univ` ⇒ `∇v ≡ 0` and `ω ≡ 0` (`global_core_forces_rigid_motion`) | 5 |

## The v2 analytic gaps are now proved

| v2 `sorry` | v3 proof |
|---|---|
| `flat_T3_killing_rigidity` | `∂ₖ∂ⱼvᵢ` is symmetric in `(k, j)` (Mathlib `ContDiffAt.isSymmSndFDerivAt`) and antisymmetric in `(j, i)` (from `S(v) = 0`), so it is `0`. Hence `v` is affine (`is_const_of_fderiv_eq_zero`), and a ℤ³-periodic affine field is constant |
| `const_of_grad_eq_zero` | Mathlib `is_const_of_fderiv_eq_zero` |
| `gronwall_scalar` | `t ↦ f(t)e^{−Kt}` is antitone (`antitoneOn_of_deriv_nonpos`) |

## What remains an explicit hypothesis

* `TZTPostulates`, including the choice of the core region `K`
* for the ceiling: `CriticalRegion Λ ⊆ K`
* `hconst` and `hinv` in `H3_mean_conserved`: a spatially constant axis and anchor invariance
  `X̄(t) = X̄(0)`. Under a localized core, `S(𝒳) = S(v)` vanishes only on `K`, so
  constancy of `𝒳` no longer follows from the postulates.
* `EnstrophyBudget`: the enstrophy identity (the flux `∇·Π` integrates out) and the Hölder bound

Theorem 5 constructs models of every hypothesis bundle, so the bundles are jointly satisfiable
and no contradiction can be derived from them.

## How to read the ceiling

`ceiling_iff` shows that, for C² periodic divergence-free flows, the localized postulate with
core `CriticalRegion Λ` holds exactly when `|ω| ≤ Λ` everywhere. `criticalRegion_empty`
shows that the critical region is then empty. So the vorticity bound is not derived from
anything weaker: it is equivalent to the postulate. The formalization verifies

    postulate (core ⊇ CriticalRegion Λ)  ⇔  |ω| ≤ Λ  ⇒  BKM integral ≤ ΛT.

It does not show that a Navier–Stokes solution satisfies the postulate. The PDE itself is not
formalized here.

## Changes from v1 (carried over from v2)

* The lakefile now loads, `ℝ` and the analysis come from Mathlib, and `λ` is renamed `lam`.
* The inconsistent v1 axioms (`T3_no_homothety`, `H3_bound`, `lambda_from_incompressibility`)
  are now proved theorems with correct hypotheses. The 14 axioms that concluded `True` are removed.
* `viscosity_swallows_half` now assumes `M ≠ 0`.
* `ω := ∇×v` is a definition, not an independent field.

## Scope

* Spatial integrals are not formalized. H2 and the ceiling are stated pointwise, which is stronger.
* Step 3 (the L² envelope) is conditional on `EnstrophyBudget`. The Hˢ/Sobolev route to BKM is
  not formalized; `BKM_sup_bound` gives the BKM bound directly from the ceiling.

---
## License & Copyright

**Copyright (c) 2026 Jung Soo Kim (Ark Project).**  
**SPDX-License-Identifier: CC-BY-NC-4.0**

The main manuscript and all enclosed supplementary materials (including analytical proofs, formal verification code, CAS scripts, DNS benchmark data, and solvers) are licensed under a **[Creative Commons Attribution-NonCommercial 4.0 International License](https://creativecommons.org/licenses/by-nc/4.0/)**.