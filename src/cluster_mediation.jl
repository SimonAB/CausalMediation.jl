# Mathematical source:
# Ham, D., Cheng, C. & Li, F. (2026). Semiparametric causal mediation analysis
# in cluster-randomized trials with post-treatment confounding.
# arXiv:2611.11903v1. https://doi.org/10.48550/arXiv.2611.11903
# https://arxiv.org/abs/2611.11903
# Targets: equations (1)--(4); canonical scores: Theorem 2, equations (5)--(7);
# weighted intercept fluctuation: Appendix I.2 (pp. 53--54).
# Independent Julia implementation, not a translation of the authors' R source.
# Implementation comparison: https://github.com/daeyoungspam-stat/CRT-interventional-mediation
# Further scope, provenance and validation: docs/src/cluster-mediation.md.

"""
    ClusterMediationSpec(treatment, outcome; cluster, mediator, confounder,
                        covariates=Symbol[], π=0.5, weighting=:cluster)

Static Bernoulli cluster-randomised mediation with one continuous mediator and
one measured post-treatment confounder per member. `π` is the known cluster
assignment probability. `weighting=:cluster` sets W=1; `:individual` sets W=N.
Baseline `covariates` must be complete finite numeric columns.

The stochastic mediator draw is independent of the outcome arm's confounder
given baseline information and cluster size. It preserves dependence *within*
the mediator vector, except in the mixed own/peer intervention defining τ.
This is a separate target from the existing `MediationSpec` nested sampler.
Identification requires no cross-cluster interference, consistency, measured
mediator-outcome confounding control and joint positivity. No graph certificate
is inferred. See [`run_cluster_mediation`](@ref).

# Reference

Ham, D., Cheng, C. & Li, F. (2026), *Semiparametric causal mediation analysis
in cluster-randomized trials with post-treatment confounding*, equations (1)--(4).
[arXiv:2611.11903](https://arxiv.org/abs/2611.11903),
[doi:10.48550/arXiv.2611.11903](https://doi.org/10.48550/arXiv.2611.11903).
"""
struct ClusterMediationSpec
    treatment::Symbol
    outcome::Symbol
    cluster::Symbol
    mediator::Symbol
    confounder::Symbol
    covariates::Vector{Symbol}
    π::Float64
    weighting::Symbol

    function ClusterMediationSpec(treatment::Symbol, outcome::Symbol;
        cluster::Symbol, mediator::Symbol, confounder::Symbol,
        covariates::Vector{Symbol} = Symbol[], π::Real = 0.5,
        weighting::Symbol = :cluster)
        isfinite(π) && 0 < π < 1 || throw(ArgumentError("π must lie strictly between zero and one"))
        weighting in (:cluster, :individual) || throw(ArgumentError("weighting must be :cluster or :individual"))
        columns = [cluster, treatment, outcome, mediator, confounder, covariates...]
        allunique(columns) || throw(ArgumentError("cluster, treatment, outcome, mediator, confounder and covariates must be distinct"))
        new(treatment, outcome, cluster, mediator, confounder, copy(covariates), Float64(π), weighting)
    end
end

"""Coefficients and member/cluster covariance components for one treatment arm."""
struct _GaussianClusterArm
    β_z::Vector{Float64}
    β_m::Vector{Float64}
    Σ_within::Matrix{Float64}
    Σ_between::Matrix{Float64}
    β_y::Vector{Float64}
    outcome_location::Vector{Float64}
    outcome_scale::Vector{Float64}
end

"""Use the untransformed polynomial basis for externally specified analytic oracles."""
_GaussianClusterArm(β_z, β_m, Σ_within, Σ_between, β_y) =
    _GaussianClusterArm(β_z, β_m, Σ_within, Σ_between, β_y, Float64[], Float64[])

"""Return 1, x, and upper-triangular quadratic monomials in mathematical order."""
function _cluster_polynomial(x)
    return [1.0; x; [x[j] * x[k] for j in eachindex(x) for k in j:length(x)]]
end

"""Integrate a quadratic polynomial using E[xⱼxₖ]=μⱼμₖ+Σⱼₖ, without sampling."""
function _cluster_quadratic_mean(β, μ, Σ)
    moments = [1.0; μ; [μ[j] * μ[k] + Σ[j, k] for j in eachindex(μ) for k in j:length(μ)]]
    return dot(β, moments)
end

"""Convert complete tabular observations to independent, sorted cluster records."""
function _cluster_records(spec, data)
    columns = [spec.cluster, spec.treatment, spec.outcome, spec.mediator, spec.confounder, spec.covariates...]
    all(c -> hasproperty(data, c), columns) || throw(ArgumentError("required cluster-mediation columns are missing"))
    nrow(data) > 0 || throw(ArgumentError("cluster mediation requires nonempty data"))
    any(ismissing, data[!, spec.cluster]) && throw(ArgumentError("cluster identifiers cannot be missing"))
    for c in columns[2:end]
        all(x -> x isa Real && isfinite(x), data[!, c]) ||
            throw(ArgumentError("column :$c must contain complete finite numeric values"))
    end
    return map(collect(groupby(data, spec.cluster; sort = true))) do block
        A = only(unique(block[!, spec.treatment])) # checked separately for an informative error
        (id = block[1, spec.cluster], A = Int(A), N = nrow(block),
            C = isempty(spec.covariates) ? zeros(nrow(block), 0) : Matrix{Float64}(block[:, spec.covariates]),
            Z = Float64.(block[!, spec.confounder]), M = Float64.(block[!, spec.mediator]),
            Y = Float64.(block[!, spec.outcome]))
    end
end

"""Linear process basis: intercept, size, member baseline and cluster-mean baseline."""
function _cluster_process_basis(O)
    C̄ = vec(mean(O.C; dims = 1))
    return hcat(ones(O.N), fill(Float64(O.N), O.N), O.C, repeat(C̄', O.N, 1))
end

"""Map [Z; M] to own/peer summaries; singleton peer summaries are zero."""
function _cluster_summary_map(N, j)
    T = zeros(4, 2N)
    T[1, j] = T[3, N + j] = 1
    if N > 1
        T[2, 1:N] .= 1 / (N - 1)
        T[4, N+1:2N] .= 1 / (N - 1)
        T[2, j] = T[4, N + j] = 0
    end
    return T
end

"""Integrate the arm-specific quadratic outcome mean over joint Z/M moments."""
function _cluster_outcome_mean(h, O, μ_z, Σ_z, μ_m, Σ_m; Σ_zm = zeros(O.N, O.N))
    μ = [μ_z; μ_m]
    Σ = [Σ_z Σ_zm; Σ_zm' Σ_m]
    C̄ = vec(mean(O.C; dims = 1))
    return map(1:O.N) do j
        T = _cluster_summary_map(O.N, j)
        baseline = [Float64(O.N); O.C[j, :]; C̄]
        μ_x = [baseline; T * μ]
        Σ_x = zeros(length(μ_x), length(μ_x))
        Σ_x[length(baseline)+1:end, length(baseline)+1:end] = T * Σ * T'
        if !isempty(h.outcome_location)
            μ_x = (μ_x - h.outcome_location) ./ h.outcome_scale
            Σ_x = Σ_x ./ (h.outcome_scale * h.outcome_scale')
        end
        _cluster_quadratic_mean(h.β_y, μ_x, Σ_x)
    end
end

"""Project a symmetric covariance component onto the specified positive cone."""
function _cluster_covariance_component(Σ; floor = 0.0)
    λ, V = eigen(Symmetric(Σ))
    return Matrix(Symmetric(V * Diagonal(max.(λ, floor)) * V'))
end

"""Fit linear marginal process means, random-intercept covariance and quadratic η."""
function _fit_cluster_arm(records)
    X = reduce(vcat, _cluster_process_basis.(records))
    Z = reduce(vcat, getproperty.(records, :Z))
    M = reduce(vcat, getproperty.(records, :M))
    β_z, β_m = pinv(X) * Z, pinv(X) * M
    Σ_member, Σ_peer = zeros(2, 2), zeros(2, 2)
    n_members, n_pairs = 0, 0
    outcome_rows = Vector{Float64}[]
    for O in records
        Xᵢ = _cluster_process_basis(O)
        E = hcat(O.Z - Xᵢ * β_z, O.M - Xᵢ * β_m)
        Σ_member .+= E' * E
        E_sum = vec(sum(E; dims = 1))
        Σ_peer .+= E_sum * E_sum' - E' * E
        n_members += O.N
        n_pairs += O.N * (O.N - 1)
        C̄ = vec(mean(O.C; dims = 1))
        for j in 1:O.N
            x = [Float64(O.N); O.C[j, :]; C̄; _cluster_summary_map(O.N, j) * [O.Z; O.M]]
            push!(outcome_rows, x)
        end
    end
    Σ_between = n_pairs == 0 ? zeros(2, 2) : _cluster_covariance_component(Σ_peer / n_pairs)
    Σ_within = _cluster_covariance_component(Σ_member / n_members - Σ_between; floor = 1e-8)
    outcome_inputs = reduce(hcat, outcome_rows)'
    outcome_location = vec(mean(outcome_inputs; dims = 1))
    outcome_scale = vec(std(outcome_inputs; dims = 1, corrected = false))
    outcome_scale[iszero.(outcome_scale)] .= 1
    design = reduce(hcat, [_cluster_polynomial((x - outcome_location) ./ outcome_scale)
        for x in eachrow(outcome_inputs)])'
    Y = reduce(vcat, getproperty.(records, :Y))
    return _GaussianClusterArm(β_z, β_m, Σ_within, Σ_between, pinv(design) * Y,
        outcome_location, outcome_scale)
end

"""Construct ϕ(Z|a,C,N), κ(M|a,C,N) and ξ(M|a,Z,C,N) from a joint Gaussian."""
function _cluster_laws(h, O)
    X = _cluster_process_basis(O)
    μ_z, μ_m = X * h.β_z, X * h.β_m
    I_N = Matrix{Float64}(I, O.N, O.N)
    J_N = ones(O.N, O.N)
    Σ_z = h.Σ_within[1, 1] * I_N + h.Σ_between[1, 1] * J_N
    Σ_m = h.Σ_within[2, 2] * I_N + h.Σ_between[2, 2] * J_N
    Σ_zm = h.Σ_within[1, 2] * I_N + h.Σ_between[1, 2] * J_N
    μ_m_z = μ_m + Σ_zm' * (Σ_z \ (O.Z - μ_z))
    Σ_m_z = Matrix(Symmetric(Σ_m - Σ_zm' * (Σ_z \ Σ_zm)))
    return (ϕ = MvNormal(μ_z, Symmetric(Σ_z)), κ = MvNormal(μ_m, Symmetric(Σ_m)),
        ξ = MvNormal(μ_m_z, Symmetric(Σ_m_z)))
end

"""Construct u₁,u₂,u₃ in Table 1 by exact integration of η against ϕ and κ."""
function _cluster_u(h, O, ϕ, κ)
    zero_covariance = zeros(O.N, O.N)
    return (
        u₁ = _cluster_outcome_mean(h, O, O.Z, zero_covariance, mean(κ), cov(κ)),
        u₂ = _cluster_outcome_mean(h, O, mean(ϕ), cov(ϕ), O.M, zero_covariance),
        u₃ = _cluster_outcome_mean(h, O, mean(ϕ), cov(ϕ), mean(κ), cov(κ)),
    )
end

"""Construct r₁,…,r₄ for independent own-M(1) and peer-M(0) draws in τ."""
function _cluster_r(h, O, ϕ, κ₁, κ₀)
    N = O.N
    r₁, r₂, r₃, r₄ = (zeros(N) for _ in 1:4)
    zero_covariance = zeros(N, N)
    for j in 1:N
        μ_m = copy(mean(κ₀))
        Σ_m = copy(cov(κ₀))
        μ_m[j] = mean(κ₁)[j]
        Σ_m[j, :] .= 0
        Σ_m[:, j] .= 0
        Σ_m[j, j] = cov(κ₁)[j, j]
        r₁[j] = _cluster_outcome_mean(h, O, O.Z, zero_covariance, μ_m, Σ_m)[j]
        r₄[j] = _cluster_outcome_mean(h, O, mean(ϕ), cov(ϕ), μ_m, Σ_m)[j]
        # r₂ retains the observed own mediator and integrates the control peers.
        μ₂, Σ₂ = copy(μ_m), copy(Σ_m)
        μ₂[j], Σ₂[j, j] = O.M[j], 0
        r₂[j] = _cluster_outcome_mean(h, O, mean(ϕ), cov(ϕ), μ₂, Σ₂)[j]
        # r₃ retains observed peers and integrates the treated own mediator.
        μ₃, Σ₃ = copy(O.M), zeros(N, N)
        μ₃[j], Σ₃[j, j] = μ_m[j], Σ_m[j, j]
        r₃[j] = _cluster_outcome_mean(h, O, mean(ϕ), cov(ϕ), μ₃, Σ₃)[j]
    end
    return (; r₁, r₂, r₃, r₄)
end

"""Theorem 2, equation (5), before W/N summation over members."""
function _cluster_θ_score(A, Y, η, ω, u, π, a, a′)
    π_a, π_a′ = a == 1 ? π : 1 - π, a′ == 1 ? π : 1 - π
    ψθ = copy(u.u₃)
    if A == a
        ψθ .+= (ω .* (Y .- η) .+ u.u₁ .- u.u₃) ./ π_a
    end
    if A == a′
        ψθ .+= (u.u₂ .- u.u₃) ./ π_a′
    end
    return ψθ
end

"""Theorem 2, equations (6)--(7), before W/N summation over members."""
function _cluster_τ_score(A, Y, η, ω, r, π)
    if A == 1
        return (ω .* (Y .- η) .+ r.r₁ .+ r.r₂ .- 2r.r₄) ./ π .+ r.r₄
    end
    return (r.r₃ .- r.r₄) ./ (1 - π) .+ r.r₄
end

"""Evaluate a density ratio on the log scale; never silently cap an uncapped ratio."""
function _cluster_ratio(log_ratio, cap)
    value = cap === nothing ? exp(log_ratio) : exp(clamp(log_ratio, -log(cap), log(cap)))
    isfinite(value) && value > 0 || throw(ArgumentError("non-finite or zero mediator density ratio; inspect joint overlap"))
    return value
end

"""Evaluate the canonical θ₀₀, θ₁₀, θ₁₁, τ scores and the ordinary trial TE score."""
function _cluster_scores(O, h, π, W, cap)
    laws = (_cluster_laws(h[1], O), _cluster_laws(h[2], O))
    η = map(h) do arm
        _cluster_outcome_mean(arm, O, O.Z, zeros(O.N, O.N), O.M, zeros(O.N, O.N))
    end
    ψ = zeros(5)
    R, H = zeros(4), zeros(4)
    log_ratios = Float64[]
    for (k, (a, a′)) in enumerate(((0, 0), (1, 0), (1, 1)))
        u = _cluster_u(h[a+1], O, laws[a+1].ϕ, laws[a′+1].κ)
        # Ratios are needed only on their residual's assigned arm.
        ω = 0.0
        if O.A == a
            log_ω = logpdf(laws[a′+1].κ, O.M) - logpdf(laws[a+1].ξ, O.M)
            push!(log_ratios, log_ω)
            ω = _cluster_ratio(log_ω, cap)
        end
        ψ[k] = W * mean(_cluster_θ_score(O.A, O.Y, η[a+1], ω, u, π, a, a′))
        if O.A == a
            π_a = a == 1 ? π : 1 - π
            H[k] = W * ω / π_a
            R[k] = H[k] * mean(O.Y - η[a+1])
        end
    end
    r = _cluster_r(h[2], O, laws[2].ϕ, laws[2].κ, laws[1].κ)
    ωτ = zeros(O.N)
    if O.A == 1
        for j in 1:O.N
            peers = [k for k in 1:O.N if k != j]
            log_κ_j = logpdf(Normal(mean(laws[2].κ)[j], sqrt(cov(laws[2].κ)[j, j])), O.M[j])
            log_κ_peers = isempty(peers) ? 0.0 : logpdf(
                MvNormal(mean(laws[1].κ)[peers], Symmetric(cov(laws[1].κ)[peers, peers])), O.M[peers])
            log_ωτ = log_κ_j + log_κ_peers - logpdf(laws[2].ξ, O.M)
            push!(log_ratios, log_ωτ)
            ωτ[j] = _cluster_ratio(log_ωτ, cap)
        end
    end
    ψ[4] = W * mean(_cluster_τ_score(O.A, O.Y, η[2], ωτ, r, π))
    if O.A == 1
        H[4] = W * mean(ωτ) / π
        R[4] = W * mean(ωτ .* (O.Y - η[2])) / π
    end
    ψ[5] = W * mean(O.Y) * (O.A / π - (1 - O.A) / (1 - π))
    n_capped = cap === nothing ? 0 : count(x -> abs(x) > log(cap), log_ratios)
    return (; ψ, R, H, min_log_ratio = minimum(log_ratios), max_log_ratio = maximum(log_ratios),
        n_capped, n_ratios = length(log_ratios))
end

"""
Solve Σᵢ(Rᵢ-Hᵢε)=0 in one validation fold and update ψᵢ by (Wᵢ-Hᵢ)ε.

Here Rᵢ=(Wᵢ/Nᵢ)ΣⱼHᵢⱼ(Yᵢⱼ-ηᵢⱼ) and Hᵢ=(Wᵢ/Nᵢ)ΣⱼHᵢⱼ include the
assignment probability. Adding ε to η and all compatible integrated means
leaves auxiliary differences unchanged. This is Appendix I.2's intercept
fluctuation, not weight truncation or a replacement uncertainty estimator.
"""
function _stabilise_cluster_fold!(ψ, R, H, W)
    denominators = vec(sum(H; dims = 1))
    all(x -> isfinite(x) && x > 0, denominators) ||
        throw(ArgumentError("fold stabilisation requires positive finite residual-weight sums"))
    ε = vec(sum(R; dims = 1)) ./ denominators
    all(isfinite, ε) || throw(ArgumentError("non-finite fold fluctuation; inspect joint overlap"))
    ψ .+= (W .- H) .* ε'
    return ε
end

"""
    run_cluster_mediation(spec, data; folds=5, ratio_cap=nothing, stabilise=false,
                          rng=Random.default_rng())

Canonical cross-fitted one-step reference estimator of Ham, Cheng & Li (2026),
Theorem 2. Linear Gaussian marginal Z/M models use member baselines, cluster
means and size; a bivariate Gaussian random intercept preserves joint mediator
dependence. Arm-specific quadratic outcome regressions use these baselines and
own/leave-one-out peer Z/M means. Gaussian integration is exact for this model.
Singleton peer means are zero. All members stay in their cluster's fold.

Return `estimates`, `se`, `covariance`, `influence`, `table`, `cluster_ids`,
`fold_ids`, and `diagnostics`. Component order is `ide, iie, iime, isme,
stochastic_total, te`; joint influence curves retain covariance of contrasts.
Diagnostics include extrema of active log density ratios before capping and
the fraction of evaluated ratios capped by the selected route. `stabilise=true`
applies Appendix I.2's weighted intercept fluctuation separately per fold,
target and weighting. `fold_fluctuations` and `fold_residual_sums` record its
coefficients and residual equations in order θ₀₀, θ₁₀, θ₁₁, τ.
`te` is a separate inverse-randomisation-weighted ordinary treatment effect.
`stochastic_total=ide+iie` need not equal `te`.

Uncapped intervals are asymptotic, require correct nuisance limits, suitable
rates, overlap and finite score moments, and are not a finite-sample coverage
guarantee. `ratio_cap≥1` opts into capped point fits and withholds *all* returned
uncertainty (`:truncated_weights`). The default does not truncate ratios.
This experimental parametric reference does not implement the paper's flexible
or classification-based estimators, missingness, longitudinal
mediation, or automatic identification. See the cluster mediation manual.

# References

Ham, D., Cheng, C. & Li, F. (2026), *Semiparametric causal mediation analysis
in cluster-randomized trials with post-treatment confounding*, Theorem 2,
equations (5)--(7), and Appendix I.2.
[arXiv:2611.11903](https://arxiv.org/abs/2611.11903),
[doi:10.48550/arXiv.2611.11903](https://doi.org/10.48550/arXiv.2611.11903).
See the manual for restrictions of this Gaussian reference implementation.
"""
function run_cluster_mediation(spec::ClusterMediationSpec, data::AbstractDataFrame;
    folds::Int = 5, ratio_cap::Union{Nothing, Real} = nothing,
    stabilise::Bool = false,
    rng::AbstractRNG = Random.default_rng())
    folds >= 2 || throw(ArgumentError("at least two cluster folds are required"))
    ratio_cap === nothing || (isfinite(ratio_cap) && ratio_cap >= 1) ||
        throw(ArgumentError("ratio_cap must be nothing or a finite number ≥ 1"))
    hasproperty(data, spec.treatment) && hasproperty(data, spec.cluster) ||
        throw(ArgumentError("treatment and cluster columns are required"))
    for block in groupby(data, spec.cluster)
        all(x -> !ismissing(x) && x in (0, 1), block[!, spec.treatment]) &&
            length(unique(block[!, spec.treatment])) == 1 ||
            throw(ArgumentError("treatment must be constant binary 0/1 within every cluster"))
    end
    records = _cluster_records(spec, data)
    K = length(records)
    fold_ids = zeros(Int, K)
    for a in 0:1
        indices = findall(O -> O.A == a, records)
        length(indices) >= folds && length(indices) - cld(length(indices), folds) >= 2 ||
            throw(ArgumentError("each treatment arm needs at least folds clusters and two training clusters per fold"))
        for (r, local_indices) in enumerate(CausalTargeted.crossfit_indices(length(indices), folds, rng))
            fold_ids[indices[local_indices]] .= r
        end
    end
    W = Float64[spec.weighting == :cluster ? 1 : O.N for O in records]
    ψ = zeros(K, 5)
    R, H = zeros(K, 4), zeros(K, 4)
    fold_fluctuations, fold_residual_sums = zeros(folds, 4), zeros(folds, 4)
    min_log_ratio, max_log_ratio = Inf, -Inf
    n_capped, n_ratios = 0, 0
    for r in 1:folds
        h = ntuple(2) do k
            _fit_cluster_arm([records[i] for i in 1:K if fold_ids[i] != r && records[i].A == k - 1])
        end
        for i in findall(==(r), fold_ids)
            scored = _cluster_scores(records[i], h, spec.π, W[i], ratio_cap)
            ψ[i, :] = scored.ψ
            R[i, :], H[i, :] = scored.R, scored.H
            min_log_ratio = min(min_log_ratio, scored.min_log_ratio)
            max_log_ratio = max(max_log_ratio, scored.max_log_ratio)
            n_capped += scored.n_capped
            n_ratios += scored.n_ratios
        end
        validation = findall(==(r), fold_ids)
        if stabilise
            fold_fluctuations[r, :] = _stabilise_cluster_fold!(
                view(ψ, validation, 1:4), view(R, validation, :), view(H, validation, :), W[validation])
        end
        fold_residual_sums[r, :] = vec(sum(
            R[validation, :] - H[validation, :] .* fold_fluctuations[r, :]'; dims = 1))
    end
    # θ̂ = Pₖψ/PₖW, D̂ᵢ = (ψᵢ - Wᵢθ̂)/PₖW.
    θ̂ = vec(sum(ψ; dims = 1)) / sum(W)
    D̂ = (ψ - W * θ̂') / mean(W)
    L = [-1.0 1 0 0 0; 0 -1 1 0 0; 0 -1 0 1 0; 0 0 1 -1 0;
        -1 0 1 0 0; 0 0 0 0 1]
    names = (:ide, :iie, :iime, :isme, :stochastic_total, :te)
    estimates = NamedTuple{names}(Tuple(L * θ̂))
    influence = D̂ * L'
    covariance = cov(influence) / K
    status = ratio_cap === nothing ? :asymptotic : :truncated_weights
    if ratio_cap !== nothing
        covariance .= NaN
        influence .= NaN
    end
    se = NamedTuple{names}(Tuple(sqrt.(max.(diag(covariance), 0))))
    critical = quantile(TDist(K - 1), 0.975)
    table = DataFrame(effect = collect(names), estimate = collect(values(estimates)),
        se = collect(values(se)), inference_status = fill(status, length(names)))
    table.lower = table.estimate - critical * table.se
    table.upper = table.estimate + critical * table.se
    diagnostics = (n_clusters = K, total_weight = sum(W), weighting = spec.weighting,
        integration = :exact_gaussian, nuisance_model = :gaussian_random_intercept_quadratic,
        inference_status = status, π = spec.π, folds = folds, ratio_cap = ratio_cap,
        min_log_ratio = min_log_ratio, max_log_ratio = max_log_ratio,
        capped_fraction = n_capped / n_ratios,
        stabilise = stabilise, fold_fluctuations = fold_fluctuations,
        fold_residual_sums = fold_residual_sums,
        identification = :user_assumed, source = "Ham, Cheng & Li (2026), Theorem 2")
    return (; spec, estimates, se, covariance, influence, table,
        cluster_ids = getproperty.(records, :id), fold_ids, diagnostics)
end
