"""
    CausalMediation

Cross-fitted mediation estimation for the CDCS stack: interventional (RI),
natural, organic, controlled direct, and recanting-twin path-specific effects.

Depends on CausalDynamics for identification certificates and CausalTargeted
for Super Learner, `ShiftPolicy`, fold helpers, and MTP density ratios.

# Documentation

- Design: `DESIGN.md`
- Boundaries: `BOUNDARIES.md`
- Naming vs R `crumble`: `NAMING.md`
"""
module CausalMediation

using DataFrames
using CausalDynamics
using CausalTargeted
using Graphs
using LinearAlgebra
using Logging
using Random
using StableRNGs
using Statistics
using StatsBase
using Distributions

include("effects.jl")
include("identify_bridge.jl")
include("nuisance.jl")
include("mediation_eif.jl")
include("eif_interventional.jl")
include("mediation_fold_cache.jl")
include("estimators.jl")
include("eif_natural.jl")
include("eif_organic.jl")
include("eif_recanting.jl")
include("eif_controlled.jl")
include("mediation_grid.jl")
include("grid.jl")
include("diagnostics.jl")
include("synthetic.jl")
include("tmle3_mediation.jl")
include("ppl_mediation.jl")
include("target_trial.jl")
include("riesz.jl")

export MediationEffect
export InterventionalMediation, NaturalMediation, OrganicMediation
export RecantingTwinMediation, ControlledDirect
export MediationSpec, MediationResult
export assumptions, assert_natural_admissible!, assert_moc_for_ri!
export plan_mediation, spec_from_identification
export run_mediation, run_mediation_grid, run_mediation_scalar
export run_mediation_scalar_ppl, run_tmle3_nde
export decompose, decompose_mediation_eif
export eif_psi_interventional, mediator_density_ratio_vs_obs
export MediationFoldCache, build_mediation_fold_cache
export prepare_ppl_mediation_spec, conjugate_mediation_bootstrap
export mediation_n_mc_sweep, mediation_stability_summary, mediation_stability_markdown
export simulate_mediation, simulate_continuous_mtp_mediation
export simulate_intermediate_confounding_mediation
export simulate_recanting_twin_mediation
export TargetTrialMediation, target_trial_mediation
export fit_riesz_representer, riesz_available

end
