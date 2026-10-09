import Lake
open Lake DSL

package «hugger_tzt» where
  version := "1.0"
  -- Ark Project Unified License (CC-BY-NC-4.0)

@[default_target]
lean_lib «HuggerTZT» where
  roots := #[
    `Theorem0_Basic_Topology,
    `Theorem1_Helmholtz_Bridge,
    `Theorem2_Integral_Alignment,
    `Theorem3_Mean_Conserved,
    `Theorem4_Viscous_Absorption
  ]