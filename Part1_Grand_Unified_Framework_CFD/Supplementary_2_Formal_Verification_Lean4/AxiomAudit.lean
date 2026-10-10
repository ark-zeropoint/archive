/-
  AxiomAudit.lean — not part of the library; run with

      lake env lean AxiomAudit.lean

  Expected output (v3): every result depends only on Lean's standard axioms
  `propext`, `Classical.choice`, `Quot.sound` — no `sorryAx`, no project axiom.
-/
import Theorem6_Vorticity_Ceiling

-- Theorem 1: Helmholtz bridge, localized core
#print axioms HUGGER.helmholtz_unique
#print axioms HUGGER.H1_of_reconstruction
#print axioms HUGGER.velocity_gradient_is_G0
#print axioms HUGGER.lambda_from_incompressibility
#print axioms HUGGER.strain_vanishes_on_core
#print axioms HUGGER.stretching_vanishes_on_core
-- Theorem 2: rigidity and H2
#print axioms HUGGER.T3_no_homothety
#print axioms HUGGER.hessian_eq_zero_of_eventually_strain_zero
#print axioms HUGGER.grad_eventually_const
#print axioms HUGGER.flat_T3_killing_rigidity
#print axioms HUGGER.remainder_vanishes_on_core_interior
#print axioms HUGGER.H2_on_core_interior
#print axioms HUGGER.vorticity_vanishes_of_global_core
-- Theorem 3: H3
#print axioms HUGGER.strain_X_eq_core
#print axioms HUGGER.strain_X_vanishes_on_core
#print axioms HUGGER.const_of_grad_eq_zero
#print axioms HUGGER.X_constant_of_killing
#print axioms HUGGER.X_constant_of_global_core
#print axioms HUGGER.H3_mean_conserved
-- Theorem 4: absorption and Grönwall
#print axioms HUGGER.viscosity_swallows_half
#print axioms HUGGER.young_enstrophy
#print axioms HUGGER.viscous_absorption
#print axioms HUGGER.gronwall_scalar
#print axioms HUGGER.enstrophy_gronwall_L2
-- Theorem 5: consistency and the flow-with-vorticity witness
#print axioms HUGGER.constState_postulates
#print axioms HUGGER.zeroBudget
#print axioms HUGGER.velocityState_postulates
#print axioms HUGGER.kolmogorov_curl
#print axioms HUGGER.kolmogorov_smoothPeriodic
#print axioms HUGGER.localized_postulates_admit_vorticity
#print axioms HUGGER.global_core_forces_rigid_motion
-- Theorem 6: vorticity ceiling
#print axioms HUGGER.vorticity_ceiling
#print axioms HUGGER.criticalRegion_empty
#print axioms HUGGER.ceiling_iff
#print axioms HUGGER.BKM_sup_bound
