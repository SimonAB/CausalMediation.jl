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
    fitted_schema.covariates == cols || throw(ArgumentError(
        "provided schema covariates $(repr(fitted_schema.covariates)) do not match $(repr(cols))",
    ))
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

"""Choose whether the implemented mediation score has a validated EIF contract."""
function _mediation_inference_status(
    a0::AbstractVector,
    a1::AbstractVector;
    estimator::Symbol,
    moc::AbstractVector{<:Symbol} = Symbol[],
    organic::Bool = false,
    ipcw_w::AbstractVector{<:Real} = ones(length(a0)),
    observed_treatment::Union{Nothing, AbstractVector} = nothing,
    binary_treatment::Bool = false,
    factor_treatment::Bool = false,
)
    length(a0) == length(a1) || throw(ArgumentError("treatment arm lengths must match"))
    length(ipcw_w) == length(a0) || throw(ArgumentError("IPCW weight length must match treatment arms"))
    all(isequal.(a0, a1)) && return :trivial
    estimator === :onestep || return :unsupported_estimator
    isempty(moc) || return :unsupported_moc
    organic && return :unsupported_estimand
    all(w -> isapprox(w, 1.0; atol = 1e-12, rtol = 0.0), ipcw_w) || return :unsupported_ipcw
    all(isequal(x, first(a0)) for x in a0) || return :unsupported_policy
    all(isequal(x, first(a1)) for x in a1) || return :unsupported_policy
    (binary_treatment || factor_treatment) || return :unsupported_policy
    if observed_treatment !== nothing
        observed = Set(observed_treatment)
        all(x -> x in observed, (first(a0), first(a1))) || return :unsupported_positivity
    end
    return :available
end

"""Replace unvalidated mediation standard errors while preserving point fits."""
function _without_mediation_inference(result, reason::Symbol)
    est, se, ic = result
    unavailable = NamedTuple{keys(se)}(ntuple(_ -> NaN, length(keys(se))))
    unavailable_ic = NamedTuple{keys(ic)}(ntuple(i -> begin
        values_i = values(ic)[i]
        values_i isa AbstractArray ? fill(NaN, size(values_i)) : NaN
    end, length(keys(ic))))
    return est, merge(unavailable, (inference_status = reason,)), unavailable_ic
end

"""Static binary-arm propensity weights with structural off-arm zeros preserved."""
function _static_binary_arm_weights(
    A::AbstractVector{<:Real},
    e::AbstractVector{<:Real};
    trunc::Real = 10.0,
)
    length(A) == length(e) || throw(ArgumentError("treatment and propensity lengths must match"))
    A_float = Float64.(A)
    e_float = Float64.(e)
    H1 = CausalTargeted.truncate_weights(A_float ./ e_float; trunc = trunc)
    H0 = CausalTargeted.truncate_weights((1 .- A_float) ./ (1 .- e_float); trunc = trunc)
    H1[A_float .== 0.0] .= 0.0
    H0[A_float .== 1.0] .= 0.0
    return H0, H1
end

"""Apply IPCW weights to mediation influence vectors when weights differ from one."""
function _summarise_mediation_influence(
    psi_nde::AbstractVector{<:Real},
    psi_nie::AbstractVector{<:Real},
    psi_te::AbstractVector{<:Real},
    ipcw_w::AbstractVector{<:Real},
    inference_status::Symbol = :available,
)
    n = length(psi_te)
    if inference_status === :trivial
        est = (nde = 0.0, nie = 0.0, te = 0.0)
        se = (nde = 0.0, nie = 0.0, te = 0.0, inference_status = inference_status)
        ic = (nde = zeros(n), nie = zeros(n), te = zeros(n))
        return est, se, ic
    end
    uniform = all(w -> isapprox(w, 1.0; atol = 1e-12, rtol = 0.0), ipcw_w)
    if uniform
        est = (nde = mean(psi_nde), nie = mean(psi_nie), te = mean(psi_te))
        se_values = (
            nde = std(psi_nde .- est.nde) / sqrt(n),
            nie = std(psi_nie .- est.nie) / sqrt(n),
            te = std(psi_te .- est.te) / sqrt(n),
        )
        ic = (nde = Float64.(psi_nde), nie = Float64.(psi_nie), te = Float64.(psi_te))
    else
        s_nde = weighted_influence_summary(psi_nde, ipcw_w)
        s_nie = weighted_influence_summary(psi_nie, ipcw_w)
        s_te = weighted_influence_summary(psi_te, ipcw_w)
        est = (nde = s_nde.estimate, nie = s_nie.estimate, te = s_te.estimate)
        se_values = (nde = s_nde.se, nie = s_nie.se, te = s_te.se)
        ic = (nde = s_nde.ic, nie = s_nie.ic, te = s_te.ic)
    end
    se = inference_status === :available ?
        merge(se_values, (inference_status = inference_status,)) :
        (nde = NaN, nie = NaN, te = NaN, inference_status = inference_status)
    if inference_status !== :available
        ic = (nde = fill(NaN, n), nie = fill(NaN, n), te = fill(NaN, n))
    end
    return est, se, ic
end
