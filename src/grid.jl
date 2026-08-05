"""Unified `run_mediation` driver and scalar API."""

using Base.Threads

"""
    run_mediation(spec, data; folds, learners, estimator, deltas, n_mc, kwargs...) -> MediationResult

Canonical mediation entry point. Dispatches on `spec.effect` and returns a
[`MediationResult`](@ref) whose `table` is the δ-grid from
[`run_mediation_grid`](@ref).

# Keyword arguments

- `estimator`: `:plugin`, `:onestep` (default), or `:tmle`
- `deltas`: MTP shift grid (defaults to CausalTargeted `default_deltas()`)
- `n_mc`: nested mediator Monte Carlo draws per unit (default `32`)
- `folds`, `learners`, `parallel`, `cache_nuisances`: Super Learner / cross-fit controls from CausalTargeted

Natural effects with nonempty `spec.moc` are refused by
[`assert_natural_admissible!`](@ref).
"""
function run_mediation(
    spec::MediationSpec,
    data::DataFrame;
    folds = mtp_settings().folds,
    learners = DEFAULT_SL_LEARNERS,
    estimator::Symbol = :onestep,
    deltas = nothing,
    n_mc::Int = 32,
    rng::AbstractRNG = StableRNG(42),
    parallel::Bool = nthreads() > 1,
    cache_nuisances::Bool = true,
    kwargs...,
)
    assert_natural_admissible!(spec)
    assert_moc_for_ri!(spec)

    δs = deltas === nothing ? default_deltas() : deltas
    table = run_mediation_grid(
        data, spec.treatment, spec.outcome;
        covar = spec.covariates,
        mediators = spec.mediators,
        moc = spec.moc,
        deltas = δs,
        folds = folds,
        learners = learners,
        estimator = estimator,
        n_mc = n_mc,
        rng = rng,
        parallel = parallel,
        cache_nuisances = cache_nuisances,
        effect = spec.effect,
        lower_q = spec.policy_d0.lower_q,
        upper_q = spec.policy_d0.upper_q,
        shift_scale = spec.policy_d0.scale,
        kwargs...,
    )
    est, se, ic = _summarise_grid(table)
    return MediationResult(
        spec, est, se, ic,
        (n_mc = n_mc, estimator = estimator, n_rows = nrow(table)),
        table,
    )
end

function _summarise_grid(table::DataFrame)
    sub = table[.!isapprox.(table.delta, 0; atol = 1e-12), :]
    isempty(sub) && (sub = table)
    d0 = first(sort(unique(Float64.(sub.delta))))
    rows = sub[Float64.(sub.delta) .== d0, :]
    get_est(lab) = begin
        r = rows[string.(rows.estimand) .== lab, :]
        isempty(r) ? (NaN, NaN) : (Float64(r.est[1]), Float64(r.se[1]))
    end
    nde, nde_se = get_est("NDE")
    nie, nie_se = get_est("NIE")
    te, te_se = get_est("TE")
    est = (nde = nde, nie = nie, te = te)
    se = (nde = nde_se, nie = nie_se, te = te_se)
    ic = (nde = Float64[], nie = Float64[], te = Float64[])
    return est, se, ic
end

function _result_table(est, se)
    rows = Dict{String, Any}[]
    for lab in keys(est)
        e = getfield(est, lab)
        s = haskey(se, lab) ? getfield(se, lab) : NaN
        lwr, upr = wald_ci(e, s)
        push!(rows, Dict(
            "effect" => uppercase(string(lab)),
            "estimate" => e,
            "se" => s,
            "lower" => lwr,
            "upper" => upr,
        ))
    end
    return DataFrame(rows)
end

"""
    run_mediation_scalar(data, trt, outcome; mediators, covar, moc, kwargs...) -> DataFrame

Binary contrast `d0=0` vs `d1=1` with NDE / NIE / TE rows.
"""
function run_mediation_scalar(
    data::DataFrame,
    trt::Symbol,
    outcome::Symbol;
    mediators::Vector{Symbol},
    covar::Vector{Symbol},
    moc::Vector{Symbol} = Symbol[],
    folds::Int = mtp_settings().folds,
    epochs::Int = 1,
    learners = DEFAULT_SL_LEARNERS,
    n_mc::Int = 32,
    estimator::Symbol = :onestep,
    effect::MediationEffect = InterventionalMediation(),
    rng = StableRNG(42),
)
    cols = unique(vcat([trt, outcome], covar, mediators, moc))
    df = dropmissing(data[:, cols])
    n = nrow(df)
    a0 = zeros(n)
    a1 = ones(n)
    est, se, _ = _effects_dispatch(
        effect, df, outcome, trt, covar, mediators, a0, a1, folds, rng;
        learners = learners,
        n_mc = n_mc,
        estimator = estimator,
        moc = moc,
        epochs = epochs,
    )
    return _result_table(est, se)
end
