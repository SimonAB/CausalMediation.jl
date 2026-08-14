"""Thin nuisance wrappers around CausalTargeted Super Learner."""

"""Fit an outcome / mediator regression Super Learner on selected columns."""
function _fit_sl_outcome(
    df::DataFrame,
    cols::Vector{Symbol},
    y::AbstractVector{<:Real};
    treatment = nothing,
    learners = DEFAULT_SL_LEARNERS,
    rng = StableRNG(1),
    schema::Union{Nothing, CausalTargeted.CovariateSchema} = nothing,
)
    fitted_schema = schema === nothing ? CausalTargeted.fit_covariate_schema(df, cols) : schema
    fitted_schema.covariates == cols || throw(ArgumentError(
        "provided schema covariates $(repr(fitted_schema.covariates)) do not match $(repr(cols))",
    ))
    X = design_matrix(fitted_schema, df; treatment = treatment)
    return fit_super_learner(X, y; learners = learners, rng = rng)
end

"""Predict from a Super Learner fit on a design matrix for `cols`."""
function _predict_sl(
    sl,
    df::DataFrame,
    cols::Vector{Symbol};
    treatment = nothing,
    treatment_values = nothing,
    schema::Union{Nothing, CausalTargeted.CovariateSchema} = nothing,
)
    fitted_schema = schema === nothing ? CausalTargeted.fit_covariate_schema(df, cols) : schema
    X = design_matrix(fitted_schema, df; treatment = treatment, treatment_values = treatment_values)
    return predict_super_learner(sl, X)
end

"""Residual SD of a mediator regression on the training fold."""
function _mediator_residual_sd(
    train::DataFrame,
    med_model,
    m_col::Symbol,
    parents,
    trt;
    schema = nothing,
)
    μ = _predict_sl(med_model, train, parents; treatment = trt, schema = schema)
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

"""Apply IPCW weights to mediation influence vectors when weights differ from one."""
function _summarise_mediation_influence(
    psi_nde::AbstractVector{<:Real},
    psi_nie::AbstractVector{<:Real},
    psi_te::AbstractVector{<:Real},
    ipcw_w::AbstractVector{<:Real},
)
    n = length(psi_te)
    uniform = all(w -> isapprox(w, 1.0; atol = 1e-12, rtol = 0.0), ipcw_w)
    if uniform
        est = (nde = mean(psi_nde), nie = mean(psi_nie), te = mean(psi_te))
        se = (
            nde = std(psi_nde .- est.nde) / sqrt(n),
            nie = std(psi_nie .- est.nie) / sqrt(n),
            te = std(psi_te .- est.te) / sqrt(n),
        )
        ic = (nde = Float64.(psi_nde), nie = Float64.(psi_nie), te = Float64.(psi_te))
        return est, se, ic
    end
    s_nde = weighted_influence_summary(psi_nde, ipcw_w)
    s_nie = weighted_influence_summary(psi_nie, ipcw_w)
    s_te = weighted_influence_summary(psi_te, ipcw_w)
    est = (nde = s_nde.estimate, nie = s_nie.estimate, te = s_te.estimate)
    se = (nde = s_nde.se, nie = s_nie.se, te = s_te.se)
    ic = (nde = s_nde.ic, nie = s_nie.ic, te = s_te.ic)
    return est, se, ic
end
