# Bayesian hierarchical Gaussian clustered stochastic mediation (Milestone 3).
# Targets follow Ham, Cheng & Li (2026), arXiv:2611.11903.
# Mechanism fits reuse CausalTargeted's structured Gaussian VI engine.

"""
    cluster_mediation_design_dimensions(spec::ClusterMediationSpec; outcome_degree=2)

Return `(confounder, mediator, outcome)` coefficient counts for explicit
`ClusterGaussianPrior` construction in [`fit_bayesian_cluster_mediation`](@ref).
"""
function cluster_mediation_design_dimensions(spec::ClusterMediationSpec; outcome_degree::Integer = 2)
    outcome_degree in (1, 2) || throw(ArgumentError("outcome_degree must be 1 (linear) or 2 (quadratic)"))
    d_base = 1 + 2 * length(spec.covariates) # [N, C_j, C̄]
    p_z = 2 + d_base                         # [1, A, N, C_j, C̄]
    p_m = 2 + d_base + 2                     # [1, A, N, C_j, C̄, Z_j, Z̄_{-j}]
    # Outcome basis: [ϕ_d(x_no_A); A .* ϕ_d(x_no_A)] where x_no_A = [N, C_j, C̄, Z_j, Z̄_{-j}, M_j, M̄_{-j}]
    d_y_in = d_base + 4
    q_y = outcome_degree == 1 ? (1 + d_y_in) : (1 + d_y_in + div(d_y_in * (d_y_in + 1), 2))
    return (; confounder = p_z, mediator = p_m, outcome = 2 * q_y)
end

function _peer_leave_one_out_mean(v::AbstractVector{Float64})
    N = length(v)
    N <= 1 && return zeros(Float64, N)
    s = sum(v)
    inv_nm1 = 1.0 / (N - 1)
    return [(s - v[j]) * inv_nm1 for j in 1:N]
end

function _bayesian_cm_basis_vector(x::AbstractVector{Float64}, degree::Int)
    return degree == 1 ? [1.0; x] : _cluster_polynomial(x)
end

function _bayesian_cm_quad_mean(β::AbstractVector{Float64}, μ::AbstractVector{Float64},
        Σ::AbstractMatrix{Float64}, degree::Int)
    if degree == 1
        return β[1] + dot(view(β, 2:length(β)), μ)
    end
    return _cluster_quadratic_mean(β, μ, Σ)
end

function _build_bayesian_cm_designs(records, outcome_degree::Int)
    rows_z = Vector{Vector{Float64}}()
    rows_m = Vector{Vector{Float64}}()
    rows_y = Vector{Vector{Float64}}()
    ys_z = Float64[]
    ys_m = Float64[]
    ys_y = Float64[]
    cids = Any[]
    for O in records
        C̄ = isempty(O.C) ? Float64[] : vec(mean(O.C; dims = 1))
        Z_peer = _peer_leave_one_out_mean(O.Z)
        M_peer = _peer_leave_one_out_mean(O.M)
        Af = Float64(O.A)
        Nf = Float64(O.N)
        for j in 1:O.N
            Cj = isempty(O.C) ? Float64[] : vec(O.C[j, :])
            base = [Nf; Cj; C̄]
            push!(rows_z, [1.0; Af; base])
            push!(rows_m, [1.0; Af; base; O.Z[j]; Z_peer[j]])
            y_in = [base; O.Z[j]; Z_peer[j]; O.M[j]; M_peer[j]]
            ϕ_y = _bayesian_cm_basis_vector(y_in, outcome_degree)
            push!(rows_y, [ϕ_y; Af .* ϕ_y])
            push!(ys_z, O.Z[j])
            push!(ys_m, O.M[j])
            push!(ys_y, O.Y[j])
            push!(cids, O.id)
        end
    end
    Xz = Matrix(reduce(hcat, rows_z)')
    Xm = Matrix(reduce(hcat, rows_m)')
    Xy = Matrix(reduce(hcat, rows_y)')
    return (; Xz, Xm, Xy, yz = ys_z, ym = ys_m, yy = ys_y, cluster_ids = cids)
end

function _fit_hierarchical_mechanism(X::Matrix{Float64}, y::Vector{Float64},
        cluster_ids::Vector, prior::CausalTargeted.ClusterGaussianPrior;
        maxiter::Int, tol::Float64, standardize::Bool)
    p = size(X, 2)
    transform = standardize ? CausalTargeted._freeze_affine_transform(X) :
        CausalTargeted.FrozenAffineDesignTransform(
            Matrix{Float64}(I, p, p), Matrix{Float64}(I, p, p),
            CausalTargeted._detect_intercept_column(X), zeros(p), ones(p))
    X_std = CausalTargeted._apply_affine_transform(X, transform)
    prior_std = CausalTargeted._transform_prior(prior, transform)
    stats_std = CausalTargeted._cluster_sufficient_stats(X_std, y, cluster_ids)
    stats_raw = CausalTargeted._cluster_sufficient_stats(X, y, cluster_ids)
    vi = CausalTargeted._structured_vi_hierarchical_gaussian(stats_std, prior_std; maxiter, tol)
    state_orig = CausalTargeted._untransform_block_state(vi.state, transform)
    return (; posterior_mean = state_orig.beta_mean,
        posterior_covariance = state_orig.beta_cov,
        block_state = state_orig,
        residual_shape = vi.residual_shape,
        residual_scale = vi.residual_scale,
        cluster_shape = vi.cluster_shape,
        cluster_scale = vi.cluster_scale,
        status = vi.status,
        converged = vi.converged,
        iterations = vi.iterations,
        elbo = vi.elbo,
        elbo_trace = vi.elbo_trace,
        prior, transform, stats_raw)
end

"""
Typed fit object returned by [`fit_bayesian_cluster_mediation`](@ref), retaining
the three hierarchical Gaussian mechanism fits (`confounder`, `mediator`,
`outcome`), `ClusterMediationSpec`, convergence status, and design metadata.
"""
struct BayesianClusterMediationFit
    spec::ClusterMediationSpec
    confounder::NamedTuple
    mediator::NamedTuple
    outcome::NamedTuple
    outcome_degree::Int
    status::Symbol
    converged::Bool
    metadata::NamedTuple
end

"""
    fit_bayesian_cluster_mediation(spec::ClusterMediationSpec, data; priors,
        outcome_degree=2, maxiter=1000, tol=1e-8, standardize=true)

Fit the factorised hierarchical Gaussian clustered mediation model
`p(Z|A,C,N) p(M|A,Z,C,N) p(Y|A,Z,M,C,N)` using `CausalTargeted`'s structured
variational engine (`q(β, b) q(σ_ε⁻²) q(σ_b⁻²)`).

`priors` must be a `NamedTuple` with explicit `ClusterGaussianPrior` entries
for `:confounder`, `:mediator`, and `:outcome` whose dimensions match
[`cluster_mediation_design_dimensions`](@ref). Separate from
[`run_cluster_mediation`](@ref); identification is marked `:user_assumed`.
"""
function fit_bayesian_cluster_mediation(spec::ClusterMediationSpec, data;
        priors::NamedTuple, outcome_degree::Integer = 2,
        maxiter::Integer = 1000, tol::Real = 1e-8, standardize::Bool = true)
    all(k -> hasproperty(priors, k), (:confounder, :mediator, :outcome)) ||
        throw(ArgumentError("priors must provide :confounder, :mediator, and :outcome ClusterGaussianPrior entries"))
    for k in (:confounder, :mediator, :outcome)
        getproperty(priors, k) isa CausalTargeted.ClusterGaussianPrior ||
            throw(ArgumentError("priors.$k must be a CausalTargeted.ClusterGaussianPrior"))
    end
    dims = cluster_mediation_design_dimensions(spec; outcome_degree)
    for k in (:confounder, :mediator, :outcome)
        length(getproperty(priors, k).mean) == getproperty(dims, k) ||
            throw(DimensionMismatch("priors.$k dimension $(length(getproperty(priors, k).mean)) must match design width $(getproperty(dims, k))"))
    end

    records = _cluster_records(spec, data)
    for O in records
        O.A in (0, 1) || throw(ArgumentError("treatment column :$(spec.treatment) must be binary 0/1"))
    end
    designs = _build_bayesian_cm_designs(records, Int(outcome_degree))

    fit_z = _fit_hierarchical_mechanism(designs.Xz, designs.yz, designs.cluster_ids,
        priors.confounder; maxiter = Int(maxiter), tol = Float64(tol), standardize)
    fit_m = _fit_hierarchical_mechanism(designs.Xm, designs.ym, designs.cluster_ids,
        priors.mediator; maxiter = Int(maxiter), tol = Float64(tol), standardize)
    fit_y = _fit_hierarchical_mechanism(designs.Xy, designs.yy, designs.cluster_ids,
        priors.outcome; maxiter = Int(maxiter), tol = Float64(tol), standardize)

    converged = fit_z.converged && fit_m.converged && fit_y.converged
    status = converged ? :converged : :maxiter_reached
    metadata = (; identification_status = :user_assumed,
        outcome_model_status = :user_assumed,
        likelihood = :hierarchical_gaussian_factorised,
        approximation = :structured_variational,
        n_clusters = length(records),
        n_observations = length(designs.yy),
        outcome_degree = Int(outcome_degree),
        dimensions = dims,
        mechanism_statuses = (; confounder = fit_z.status, mediator = fit_m.status, outcome = fit_y.status))
    return BayesianClusterMediationFit(spec, fit_z, fit_m, fit_y, Int(outcome_degree),
        status, converged, metadata)
end

function _sample_mechanism_global(mech, n_draws::Int, rng::AbstractRNG)
    p = length(mech.posterior_mean)
    Z = randn(rng, p, n_draws)
    ldiv!(mech.block_state.beta_precision_chol.U, Z)
    beta_draws = Z .+ mech.posterior_mean
    # σ² ~ InvGamma(shape, scale) ≡ 1 / Gamma(shape, 1 / scale)
    σ2_eps = [rand(rng, InverseGamma(mech.residual_shape, mech.residual_scale)) for _ in 1:n_draws]
    σ2_b = [rand(rng, InverseGamma(mech.cluster_shape, mech.cluster_scale)) for _ in 1:n_draws]
    return (; beta_draws, σ2_eps, σ2_b)
end

function _bayesian_cm_arm_beta(β_y::AbstractVector{Float64}, a_y::Float64)
    q = div(length(β_y), 2)
    return a_y == 0.0 ? view(β_y, 1:q) : (view(β_y, 1:q) .+ view(β_y, (q + 1):(2 * q)))
end

"""
Compute exact cluster-level outcome expectation for one draw of parameters
given `(μ_z, Σ_z, μ_m, Σ_m, Σ_zm)` under outcome treatment arm `a_y`.
"""
function _bayesian_cm_cluster_mean(β_y::AbstractVector{Float64}, degree::Int,
        a_y::Float64, N::Int, C::Matrix{Float64}, C̄::Vector{Float64},
        μ_z::Vector{Float64}, Σ_z::Matrix{Float64},
        μ_m::Vector{Float64}, Σ_m::Matrix{Float64}, Σ_zm::Matrix{Float64})
    β_arm = _bayesian_cm_arm_beta(β_y, a_y)
    μ = [μ_z; μ_m]
    Σ = [Σ_z Σ_zm; Σ_zm' Σ_m]
    Nf = Float64(N)
    total = 0.0
    for j in 1:N
        T = _cluster_summary_map(N, j)
        Cj = isempty(C) ? Float64[] : vec(C[j, :])
        prefix = [Nf; Cj; C̄]
        np = length(prefix)
        μ_x = [prefix; T * μ]
        Σ_x = zeros(Float64, np + 4, np + 4)
        Σ_x[(np + 1):end, (np + 1):end] = T * Σ * T'
        total += _bayesian_cm_quad_mean(β_arm, μ_x, Σ_x, degree)
    end
    return total / N
end

"""
Exact cluster-level mediation functionals for one joint parameter draw.
Returns `(te, stochastic_total, ide, iie, iime, isme, θ00, θ10, θ11, τ, μ_nat0, μ_nat1)`.
"""
function _evaluate_bayesian_cm_draw(target_clusters, degree::Int, weighting::Symbol,
        β_z::AbstractVector{Float64}, σ2_ez::Float64, σ2_bz::Float64,
        β_m::AbstractVector{Float64}, σ2_em::Float64, σ2_bm::Float64,
        β_y::AbstractVector{Float64})
    K = length(target_clusters)
    weights = weighting === :cluster ? fill(1.0 / K, K) :
        begin
            N_tot = sum(c -> c.N, target_clusters)
            [c.N / N_tot for c in target_clusters]
        end

    γ_own = β_m[end - 1]
    γ_peer = β_m[end]
    β_m_base = view(β_m, 1:(length(β_m) - 2))

    θ00 = 0.0; θ10 = 0.0; θ11 = 0.0; τ_val = 0.0
    μ_nat0 = 0.0; μ_nat1 = 0.0

    for (idx, cl) in enumerate(target_clusters)
        w = weights[idx]
        N = cl.N
        I_N = Matrix{Float64}(I, N, N)
        J_N = ones(Float64, N, N)
        zero_N = zeros(Float64, N, N)

        # Confounder laws Z | a
        μ_z0 = cl.X_base0 * β_z
        μ_z1 = cl.X_base1 * β_z
        Σ_z = σ2_ez .* I_N .+ σ2_bz .* J_N

        # Linear operator G_M mapping Z -> γ_own Z_j + γ_peer Z̄_{-j}
        G_M = N == 1 ? (γ_own .* I_N) :
            ((γ_own - γ_peer / (N - 1)) .* I_N .+ (γ_peer / (N - 1)) .* J_N)

        # Marginal mediator laws M | a (integrating out Z | a)
        μ_m0 = cl.X_base0 * β_m_base + G_M * μ_z0
        μ_m1 = cl.X_base1 * β_m_base + G_M * μ_z1
        Σ_m_given_z = σ2_em .* I_N .+ σ2_bm .* J_N
        Σ_zm = Σ_z * G_M'
        Σ_m = Matrix(Symmetric(Σ_m_given_z + G_M * Σ_z * G_M'))

        # θ(0, 0), θ(1, 0), θ(1, 1) with independent Z and M draws (Σ_zm = 0)
        θ00 += w * _bayesian_cm_cluster_mean(β_y, degree, 0.0, N, cl.C, cl.C̄, μ_z0, Σ_z, μ_m0, Σ_m, zero_N)
        θ10 += w * _bayesian_cm_cluster_mean(β_y, degree, 1.0, N, cl.C, cl.C̄, μ_z1, Σ_z, μ_m0, Σ_m, zero_N)
        θ11 += w * _bayesian_cm_cluster_mean(β_y, degree, 1.0, N, cl.C, cl.C̄, μ_z1, Σ_z, μ_m1, Σ_m, zero_N)

        # Natural joint means (retaining Σ_zm)
        μ_nat0 += w * _bayesian_cm_cluster_mean(β_y, degree, 0.0, N, cl.C, cl.C̄, μ_z0, Σ_z, μ_m0, Σ_m, Σ_zm)
        μ_nat1 += w * _bayesian_cm_cluster_mean(β_y, degree, 1.0, N, cl.C, cl.C̄, μ_z1, Σ_z, μ_m1, Σ_m, Σ_zm)

        # τ: outcome arm 1, Z ~ p(Z|1), independent own M_j ~ p(M_j|1) and peers M_{-j} ~ p(M_{-j}|0)
        β_arm1 = _bayesian_cm_arm_beta(β_y, 1.0)
        τ_cluster = 0.0
        Nf = Float64(N)
        for j in 1:N
            μ_m_tau = copy(μ_m0)
            Σ_m_tau = copy(Σ_m)
            μ_m_tau[j] = μ_m1[j]
            Σ_m_tau[j, :] .= 0.0
            Σ_m_tau[:, j] .= 0.0
            Σ_m_tau[j, j] = Σ_m[j, j]
            T = _cluster_summary_map(N, j)
            Cj = isempty(cl.C) ? Float64[] : vec(cl.C[j, :])
            prefix = [Nf; Cj; cl.C̄]
            np = length(prefix)
            μ_x = [prefix; T * [μ_z1; μ_m_tau]]
            Σ_x = zeros(Float64, np + 4, np + 4)
            Σ_joint = [Σ_z zero_N; zero_N Σ_m_tau]
            Σ_x[(np + 1):end, (np + 1):end] = T * Σ_joint * T'
            τ_cluster += _bayesian_cm_quad_mean(β_arm1, μ_x, Σ_x, degree)
        end
        τ_val += w * (τ_cluster / N)
    end

    ide = θ10 - θ00
    iie = θ11 - θ10
    stochastic_total = θ11 - θ00
    iime = τ_val - θ10
    isme = θ11 - τ_val
    te = μ_nat1 - μ_nat0
    return (; te, stochastic_total, ide, iie, iime, isme, θ00, θ10, θ11, τ = τ_val, μ_nat0, μ_nat1)
end

function _target_cluster_records(spec::ClusterMediationSpec, target)
    cols = [spec.cluster, spec.covariates...]
    all(c -> hasproperty(target, c), cols) ||
        throw(ArgumentError("target data missing required cluster/covariate columns"))
    nrow(target) > 0 || throw(ArgumentError("target data must be nonempty"))
    any(ismissing, target[!, spec.cluster]) && throw(ArgumentError("target cluster column must not contain missing values"))
    for c in spec.covariates
        all(x -> x isa Real && isfinite(x), target[!, c]) ||
            throw(ArgumentError("target covariate :$c must contain complete finite numeric values"))
    end
    return map(collect(groupby(target, spec.cluster; sort = true))) do block
        N = nrow(block)
        C = isempty(spec.covariates) ? zeros(Float64, N, 0) : Matrix{Float64}(block[:, spec.covariates])
        C̄ = isempty(spec.covariates) ? Float64[] : vec(mean(C; dims = 1))
        Nf = Float64(N)
        X_base0 = Matrix(reduce(hcat, [[1.0; 0.0; Nf; (isempty(C) ? Float64[] : vec(C[j, :])); C̄] for j in 1:N])')
        X_base1 = Matrix(reduce(hcat, [[1.0; 1.0; Nf; (isempty(C) ? Float64[] : vec(C[j, :])); C̄] for j in 1:N])')
        (; id = block[1, spec.cluster], N, C, C̄, X_base0, X_base1)
    end
end

"""
    bayesian_mediation_effects(fit::BayesianClusterMediationFit, target;
        n_draws=2000, rng=Random.default_rng(), level=0.95,
        weighting=fit.spec.weighting, correction=:none)

Compute joint posterior summaries for `te`, `stochastic_total`, `ide`, `iie`,
`iime`, and `isme` from a [`BayesianClusterMediationFit`](@ref), using shared
posterior draws across all contrasts and exact Gaussian integration over both
cluster random effects and `(Z, M)`.

If `fit.status !== :converged`, intervals and joint covariance are withheld
(`NaN` / `intervals_available = false`).
"""
function bayesian_mediation_effects(fit::BayesianClusterMediationFit, target;
        n_draws::Integer = 2000, rng::AbstractRNG = Random.default_rng(),
        level::Real = 0.95, weighting::Symbol = fit.spec.weighting,
        correction::Symbol = :none)
    correction === :none || throw(ArgumentError(
        "correction=:$correction is not yet supported on BayesianClusterMediationFit; use correction=:none"))
    weighting in (:cluster, :individual) || throw(ArgumentError("weighting must be :cluster or :individual"))
    n_draws >= 2 || throw(ArgumentError("n_draws must be at least 2"))
    isfinite(level) && 0 < level < 1 || throw(ArgumentError("level must lie strictly between zero and one"))

    target_clusters = _target_cluster_records(fit.spec, target)
    effect_names = [:te, :stochastic_total, :ide, :iie, :iime, :isme]
    m = length(effect_names)
    S = Int(n_draws)

    if fit.status !== :converged
        # Evaluate point estimate at variational posterior means only; withhold uncertainty
        pt = _evaluate_bayesian_cm_draw(target_clusters, fit.outcome_degree, weighting,
            fit.confounder.posterior_mean,
            fit.confounder.residual_scale / (fit.confounder.residual_shape - 1.0),
            fit.confounder.cluster_scale / (fit.confounder.cluster_shape - 1.0),
            fit.mediator.posterior_mean,
            fit.mediator.residual_scale / (fit.mediator.residual_shape - 1.0),
            fit.mediator.cluster_scale / (fit.mediator.cluster_shape - 1.0),
            fit.outcome.posterior_mean)
        est = [getproperty(pt, e) for e in effect_names]
        nan_vec = fill(NaN, m)
        estimates = DataFrame(effect = effect_names, mean = est, sd = nan_vec,
            lower = nan_vec, upper = nan_vec, mcse = nan_vec)
        metadata = (; fit.metadata..., weighting, n_target_clusters = length(target_clusters),
            level = Float64(level), n_draws = S, correction,
            uncertainty = :unavailable, intervals = :unavailable, intervals_available = false)
        return (; estimates, covariance = fill(NaN, m, m), draws = fill(NaN, m, 0), metadata)
    end

    dz = _sample_mechanism_global(fit.confounder, S, rng)
    dm = _sample_mechanism_global(fit.mediator, S, rng)
    dy = _sample_mechanism_global(fit.outcome, S, rng)

    draws = Matrix{Float64}(undef, m, S)
    for s in 1:S
        ev = _evaluate_bayesian_cm_draw(target_clusters, fit.outcome_degree, weighting,
            view(dz.beta_draws, :, s), dz.σ2_eps[s], dz.σ2_b[s],
            view(dm.beta_draws, :, s), dm.σ2_eps[s], dm.σ2_b[s],
            view(dy.beta_draws, :, s))
        draws[1, s] = ev.te
        draws[2, s] = ev.stochastic_total
        draws[3, s] = ev.ide
        draws[4, s] = ev.iie
        draws[5, s] = ev.iime
        draws[6, s] = ev.isme
    end

    α = (1.0 - Float64(level)) / 2.0
    means = vec(mean(draws; dims = 2))
    sds = [std(view(draws, j, :)) for j in 1:m]
    mcses = sds ./ sqrt(S)
    lowers = [quantile(view(draws, j, :), α) for j in 1:m]
    uppers = [quantile(view(draws, j, :), 1.0 - α) for j in 1:m]
    cov_mat = Matrix(Symmetric(cov(draws; dims = 2)))

    estimates = DataFrame(effect = effect_names, mean = means, sd = sds,
        lower = lowers, upper = uppers, mcse = mcses)
    metadata = (; fit.metadata..., weighting, n_target_clusters = length(target_clusters),
        level = Float64(level), n_draws = S, correction,
        uncertainty = :structured_variational, intervals = :pointwise, intervals_available = true)
    return (; estimates, covariance = cov_mat, draws, metadata)
end
