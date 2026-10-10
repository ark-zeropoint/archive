import Lake
open Lake DSL

/-
  H.U.G.G.E.R + TZT — Lean 4 formal skeleton (Supplementary 2)
  Ark Project Unified License (CC-BY-NC-4.0)
-/
package «hugger_tzt» where
  leanOptions := #[⟨`autoImplicit, false⟩]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.34.1"

@[default_target]
lean_lib «HuggerTZT» where
  roots := #[
    `Theorem0_Basic_Topology,
    `Theorem1_Helmholtz_Bridge,
    `Theorem2_Integral_Alignment,
    `Theorem3_Mean_Conserved,
    `Theorem4_Viscous_Absorption,
    `Theorem5_Consistency_Audit,
    `Theorem6_Vorticity_Ceiling
  ]
