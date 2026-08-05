"""Interventional EIF helpers (moc-aware nested expectations)."""

"""
    _nested_mediator_outcome_means(...) -> (y_a0_m0, y_a1_m0, y_a1_m1)

Average outcome predictions over nested Gaussian draws of mediators (and
optional intermediate confounders) under intervened treatment values.
"""
function _nested_mediator_outcome_means(
    ols_y,
    block::DataFrame,
    adjust::Vector{Symbol},
    mediators::Vector{Symbol},
    med_models,
    σ_m::Vector{Float64},
    med_parents::Vector{Symbol},
    trt::Symbol,
    a0::AbstractVector{<:Real},
    a1::AbstractVector{<:Real},
    n_mc::Int,
    rng::AbstractRNG;
    moc::Vector{Symbol} = Symbol[],
    moc_models = nothing,
    σ_z::Vector{Float64} = Float64[],
    moc_parents::Vector{Symbol} = Symbol[],
)
    n_te = nrow(block)
    n_med = length(mediators)
    y_a0_m0 = zeros(n_te)
    y_a1_m0 = zeros(n_te)
    y_a1_m1 = zeros(n_te)

    has_moc = !isempty(moc) && moc_models !== nothing

    # Mediator means at a0 / a1 (moc columns in block used as-is for conditioning
    # when present; nested draws overwrite moc under each policy when has_moc).
    function _μ_med(a_pol)
        return hcat([
            _predict_sl(mm, block, med_parents; treatment = trt, treatment_values = a_pol)
            for mm in med_models
        ]...)
    end

    if !has_moc
        μ0 = _μ_med(a0)
        μ1 = _μ_med(a1)
        if n_mc <= 1
            block_m0 = copy(block)
            block_m1 = copy(block)
            for j in 1:n_med
                block_m0[!, mediators[j]] = μ0[:, j]
                block_m1[!, mediators[j]] = μ1[:, j]
            end
            y_a0_m0 .= _predict_sl(ols_y, block_m0, adjust; treatment = trt, treatment_values = a0)
            y_a1_m0 .= _predict_sl(ols_y, block_m0, adjust; treatment = trt, treatment_values = a1)
            y_a1_m1 .= _predict_sl(ols_y, block_m1, adjust; treatment = trt, treatment_values = a1)
            return y_a0_m0, y_a1_m0, y_a1_m1
        end
        block_m0 = copy(block)
        block_m1 = copy(block)
        block_m0_anti = copy(block)
        block_m1_anti = copy(block)
        for _ in 1:n_mc
            for j in 1:n_med
                noise0 = σ_m[j] .* randn(rng, n_te)
                noise1 = σ_m[j] .* randn(rng, n_te)
                block_m0[!, mediators[j]] = μ0[:, j] .+ noise0
                block_m1[!, mediators[j]] = μ1[:, j] .+ noise1
                block_m0_anti[!, mediators[j]] = μ0[:, j] .- noise0
                block_m1_anti[!, mediators[j]] = μ1[:, j] .- noise1
            end
            y_a0_m0 .+= _predict_sl(ols_y, block_m0, adjust; treatment = trt, treatment_values = a0) .+
                         _predict_sl(ols_y, block_m0_anti, adjust; treatment = trt, treatment_values = a0)
            y_a1_m0 .+= _predict_sl(ols_y, block_m0, adjust; treatment = trt, treatment_values = a1) .+
                         _predict_sl(ols_y, block_m0_anti, adjust; treatment = trt, treatment_values = a1)
            y_a1_m1 .+= _predict_sl(ols_y, block_m1, adjust; treatment = trt, treatment_values = a1) .+
                         _predict_sl(ols_y, block_m1_anti, adjust; treatment = trt, treatment_values = a1)
        end
        inv_mc = 1 / (2 * n_mc)
        y_a0_m0 .*= inv_mc
        y_a1_m0 .*= inv_mc
        y_a1_m1 .*= inv_mc
        return y_a0_m0, y_a1_m0, y_a1_m1
    end

    # With moc: draw Z under a_m then M | A=a_m, Z, W; evaluate Q(a_t, M, Z, W)
    n_z = length(moc)
    block_work = copy(block)
    block_anti = copy(block)
    n_rep = max(n_mc, 1)
    for _ in 1:n_rep
        for (a_m, a_t, acc) in ((a0, a0, y_a0_m0), (a0, a1, y_a1_m0), (a1, a1, y_a1_m1))
            μ_z = hcat([
                _predict_sl(zm, block, moc_parents; treatment = trt, treatment_values = a_m)
                for zm in moc_models
            ]...)
            for j in 1:n_z
                noise = σ_z[j] .* randn(rng, n_te)
                block_work[!, moc[j]] = μ_z[:, j] .+ noise
                block_anti[!, moc[j]] = μ_z[:, j] .- noise
            end
            μ_m = hcat([
                _predict_sl(mm, block_work, med_parents; treatment = trt, treatment_values = a_m)
                for mm in med_models
            ]...)
            μ_m_anti = hcat([
                _predict_sl(mm, block_anti, med_parents; treatment = trt, treatment_values = a_m)
                for mm in med_models
            ]...)
            for j in 1:n_med
                nm = σ_m[j] .* randn(rng, n_te)
                block_work[!, mediators[j]] = μ_m[:, j] .+ nm
                block_anti[!, mediators[j]] = μ_m_anti[:, j] .- nm
            end
            acc .+= _predict_sl(ols_y, block_work, adjust; treatment = trt, treatment_values = a_t) .+
                    _predict_sl(ols_y, block_anti, adjust; treatment = trt, treatment_values = a_t)
        end
    end
    inv_mc = 1 / (2 * n_rep)
    y_a0_m0 .*= inv_mc
    y_a1_m0 .*= inv_mc
    y_a1_m1 .*= inv_mc
    return y_a0_m0, y_a1_m0, y_a1_m1
end
