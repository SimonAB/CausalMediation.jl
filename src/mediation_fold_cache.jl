"""Fold-level nuisance cache for mediation δ-grids (outcome, mediators, moc, exposure)."""

"""
    MediationFoldCache

Per-fold SuperLearner fits reused across δ values. Policy-specific predictions
still depend on `a_nat` / `a_shift`.
"""
struct MediationFoldCache
    fold_sets::Vector{Vector{Int}}
    outcome_models::Vector{SuperLearnerFit}
    mediator_models::Vector{Vector{SuperLearnerFit}}
    sigma_m::Vector{Vector{Float64}}
    moc_models::Vector{Union{Nothing, Vector{SuperLearnerFit}}}
    sigma_z::Vector{Vector{Float64}}
    exposure_models::Vector{SuperLearnerFit}
    adjust::Vector{Symbol}
    covar::Vector{Symbol}
    mediators::Vector{Symbol}
    moc::Vector{Symbol}
    trt::Symbol
    learners::Tuple
    rng_seed::UInt
end

"""
    build_mediation_fold_cache(df, outcome, trt, covar, mediators, folds, rng; learners, moc) -> MediationFoldCache

Cross-fit Super Learner fits for outcome, mediators, optional `moc`, and
exposure, reused across δ values on a mediation grid.
"""
function build_mediation_fold_cache(
    df::DataFrame,
    outcome::Symbol,
    trt::Symbol,
    covar::Vector{Symbol},
    mediators::Vector{Symbol},
    folds::Int,
    rng::AbstractRNG;
    learners = DEFAULT_SL_LEARNERS,
    moc::Vector{Symbol} = Symbol[],
)
    n = nrow(df)
    y = Float64.(df[!, outcome])
    a = Float64.(df[!, trt])
    med_parents = _mediator_parents(covar, moc)
    moc_parents = copy(covar)
    adjust = _outcome_parents(covar, moc, mediators)
    fold_sets = crossfit_indices(n, folds, rng)
    seed = UInt(mod(hash(rng), typemax(UInt)))

    outcome_models = SuperLearnerFit[]
    mediator_models = Vector{SuperLearnerFit}[]
    sigma_m = Vector{Float64}[]
    moc_models_v = Union{Nothing, Vector{SuperLearnerFit}}[]
    sigma_z = Vector{Float64}[]
    exposure_models = SuperLearnerFit[]

    for test_idx in fold_sets
        train_idx = setdiff(1:n, test_idx)
        train = df[train_idx, :]
        y_tr = y[train_idx]
        ols_y = _fit_sl_outcome(train, adjust, y_tr; treatment = trt, learners = learners, rng = rng)
        med_models = SuperLearnerFit[
            _fit_sl_outcome(train, med_parents, Float64.(train[!, m]); treatment = trt, learners = learners, rng = rng)
            for m in mediators
        ]
        σ = [
            _mediator_residual_sd(train, med_models[j], mediators[j], med_parents, trt)
            for j in eachindex(mediators)
        ]
        if !isempty(moc)
            zm = SuperLearnerFit[
                _fit_sl_outcome(train, moc_parents, Float64.(train[!, z]); treatment = trt, learners = learners, rng = rng)
                for z in moc
            ]
            σz = [
                _mediator_residual_sd(train, zm[j], moc[j], moc_parents, trt)
                for j in eachindex(moc)
            ]
            push!(moc_models_v, zm)
            push!(sigma_z, σz)
        else
            push!(moc_models_v, nothing)
            push!(sigma_z, Float64[])
        end
        sl_a = fit_super_learner(
            design_matrix(train, covar), a[train_idx];
            learners = learners, rng = rng,
        )
        push!(outcome_models, ols_y)
        push!(mediator_models, med_models)
        push!(sigma_m, σ)
        push!(exposure_models, sl_a)
    end

    return MediationFoldCache(
        fold_sets, outcome_models, mediator_models, sigma_m, moc_models_v, sigma_z,
        exposure_models, adjust, covar, mediators, moc, trt, Tuple(learners), seed,
    )
end
