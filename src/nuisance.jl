"""Thin nuisance wrappers around CausalTargeted Super Learner."""

"""Fit an outcome / mediator regression Super Learner on selected columns."""
function _fit_sl_outcome(
    df::DataFrame,
    cols::Vector{Symbol},
    y::AbstractVector{<:Real};
    treatment = nothing,
    learners = DEFAULT_SL_LEARNERS,
    rng = StableRNG(1),
)
    X = design_matrix(df, cols; treatment = treatment)
    return fit_super_learner(X, y; learners = learners, rng = rng)
end

"""Predict from a Super Learner fit on a design matrix for `cols`."""
function _predict_sl(sl, df::DataFrame, cols::Vector{Symbol}; treatment = nothing, treatment_values = nothing)
    X = design_matrix(df, cols; treatment = treatment, treatment_values = treatment_values)
    return predict_super_learner(sl, X)
end

"""Residual SD of a mediator regression on the training fold."""
function _mediator_residual_sd(train::DataFrame, med_model, m_col::Symbol, parents, trt)
    μ = _predict_sl(med_model, train, parents; treatment = trt)
    return robust_residual_sd(Float64.(train[!, m_col]) .- μ)
end

"""Parents for mediator density: baseline + moc (intermediate confounders)."""
function _mediator_parents(covar::Vector{Symbol}, moc::Vector{Symbol})
    return unique(vcat(covar, moc))
end

"""Parents for outcome regression: baseline + moc + mediators."""
function _outcome_parents(covar::Vector{Symbol}, moc::Vector{Symbol}, mediators::Vector{Symbol})
    return unique(vcat(covar, moc, mediators))
end
