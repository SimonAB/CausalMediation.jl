"""tmle3-style interventional direct effect (NDE) for binary treatment."""

"""
    run_tmle3_nde(data, treatment, outcome; baseline, mediators, folds, rng) -> DataFrame

Cross-fitted SuperLearner interventional NDE with one-step TMLE targeting.
"""
function run_tmle3_nde(
    data::DataFrame,
    treatment::Symbol,
    outcome::Symbol;
    baseline::Vector{Symbol},
    mediators::Vector{Symbol} = Symbol[],
    folds::Int = mtp_settings().folds,
    rng = StableRNG(42),
    handle_missing::Symbol = :drop,
)
    all_cols = unique(vcat(baseline, mediators, [treatment]))
    df, ipcw_w, extra_cols = handle_missing_data(data, outcome, all_cols, handle_missing; rng = rng)
    baseline = columns_present(df, unique(vcat(baseline, extra_cols)))
    mediators = columns_present(df, mediators)
    n = nrow(df)
    A = Float64.(df[!, treatment])
    Y = Float64.(df[!, outcome])
    covar = columns_present(df, baseline)

    e_hat = crossfit_propensity(df, treatment, covar, folds, rng)
    adjust = isempty(mediators) ? covar : unique(vcat(covar, mediators))

    Q1 = crossfit_predict_outcome(
        df, outcome, treatment, adjust, ones(n), folds, rng,
    )
    Q0 = crossfit_predict_outcome(
        df, outcome, treatment, adjust, zeros(n), folds, rng,
    )
    Q_obs = crossfit_outcome_predictions(df, outcome, treatment, adjust, folds, rng)

    H1 = A ./ e_hat
    H0 = (1 .- A) ./ (1 .- e_hat)
    resid = Y .- Q_obs
    denom = sum((H1 .- H0) .^ 2)
    ε = denom > 1e-12 ? sum((H1 .- H0) .* resid) / denom : 0.0
    ic = (Q1 .- Q0) .+ ε .* (H1 .- H0)
    if _uses_ipcw_weights(ipcw_w)
        s = _weighted_influence_summary(ic, ipcw_w)
        est, se = s.estimate, s.se
    else
        est = mean(ic)
        se = std(ic) / sqrt(n)
    end

    return DataFrame(
        treatment = [string(treatment)],
        outcome = [string(outcome)],
        estimate = [est],
        se = [se],
    )
end
