# Clustered mediation

`ClusterMediationSpec` and `run_cluster_mediation` provide an experimental
frequentist parametric reference for the canonical Theorem 2 EIF estimator of
Ham, Cheng and Li (2026). This page also documents the separate experimental
parametric hierarchical Bayesian route (`fit_bayesian_cluster_mediation` and
`bayesian_mediation_effects`), which reuses `ClusterMediationSpec` for the same
targets. Neither route is part of the ordinary `MediationSpec` / `run_mediation`
API.

Daeyoung Ham, Chao Cheng and Fan Li (2026). *Semiparametric causal mediation
analysis in cluster-randomized trials with post-treatment confounding*.
arXiv:2611.11903v1, 9 November 2026.
[Paper](https://arxiv.org/abs/2611.11903) ·
[doi:10.48550/arXiv.2611.11903](https://doi.org/10.48550/arXiv.2611.11903).
The shared book bibliography key is `ham2026cluster`.

Equations (1)--(4) define the targets, Table 1 defines the auxiliary means,
Theorem 2 and equations (5)--(7) define the canonical scores, and Appendix I.2
(pp. 53--54) describes the optional weighted intercept fluctuation below.

This implementation is written independently from the published equations;
the authors' R source is neither translated nor distributed. Their repository's
public licensing terms were pending when this implementation was prepared.

## Target and assumptions

Treatment ``A`` is assigned to independent clusters by Bernoulli randomisation
with known common probability ``π``. Each member has a measured post-treatment
confounder ``Z_j``, continuous mediator ``M_j`` and outcome ``Y_j``. Baseline
information ``C`` includes the supplied covariates; ``N`` is cluster size.
Outcomes may depend on other members' post-treatment variables, with no
interference across clusters.

Let ``G(a′)`` be a joint mediator draw under arm ``a′``, conditional on ``C,N``,
independent of the confounder vector ``Z(a)`` retained under the outcome arm:

```math
\theta(a,a′)=\frac{E\left[\frac{W}{N}\sum_{j=1}^N
\iint\eta_j(a,z,m,C,N)\phi(z\mid a,C,N)
\kappa(m\mid a′,C,N)\,dz\,dm\right]}{E[W]}.
```

For ``τ``, replace ``κ(m|a′,C,N)`` by
``κ_j(m_j|1,C,N)κ_{-j}(m_{-j}|0,C,N)`` separately for each member ``j``.
This removes own/peer mediator dependence. Individual and spillover contrasts
therefore include dependence changes as well as marginal distribution changes;
they are not unique allocations to causal pathways. The member-specific mixed
interventions need not form one implementable assignment to the entire cluster.

| Result | Definition |
|---|---|
| `ide` | ``θ(1,0)-θ(0,0)`` |
| `iie` | ``θ(1,1)-θ(1,0)`` |
| `iime` | ``τ-θ(1,0)`` |
| `isme` | ``θ(1,1)-τ`` |
| `stochastic_total` | ``θ(1,1)-θ(0,0)=ide+iie`` |
| `te` | Ordinary treatment effect, without intervening on mediators |

`te` uses a separate inverse-randomisation-weighted score. It need not equal
`stochastic_total`. `weighting=:cluster` sets ``W=1``; `:individual` sets ``W=N``.
Cluster membership and size must be determined before treatment.

Identification requires consistency, adequate measured control of joint
mediator-outcome confounding, and joint positivity under the specified
stochastic interventions. Randomisation alone does not establish mediator-outcome
exchangeability. The runner records `identification=:user_assumed` and does not
consume the existing no-interference `MediationQuery` certificate.

## Model and equations in code

Within each training arm, marginal means for ``Z,M`` are linear in cluster size,
member baselines and cluster-mean baselines. Joint residual covariance is

```math
\operatorname{Cov}\begin{pmatrix}Z\\M\end{pmatrix}
=\Sigma_{\mathrm{within}}\otimes I_N+
\Sigma_{\mathrm{between}}\otimes\mathbf1_N\mathbf1_N^T.
```

The two 2-by-2 components are estimated from member and distinct-peer residual
cross-products and projected onto the positive semidefinite cone. Within-component
eigenvalues are bounded below by ``10^{-8}`` for numerical evaluation. This
Gaussian random-intercept reference excludes negative intracluster correlations
and arbitrary dependence. Conditional Gaussian algebra yields ``ξ(M|A,Z,C,N)``;
marginalisation yields ``ϕ(Z|A,C,N)`` and ``κ(M|A,C,N)``.

Arm-specific outcome regression uses the full quadratic polynomial in size,
member and mean baselines, and own/leave-one-out peer ``Z,M`` summaries.
Before forming quadratic monomials, outcome inputs are centred and scaled
using training-fold moments. Their Gaussian means and covariances undergo the
same affine transformation during integration. This preserves the quadratic
function space while preventing a change of confounder origin from causing
numerical rank loss. Least squares uses a pseudoinverse for structural duplicate
columns. Baseline covariates must be complete, numeric and suitably scaled for
the linear process models. Singleton peer summaries are zero.

`_cluster_u` and `_cluster_r` implement Table 1's integrated means.
`_cluster_quadratic_mean` uses ``E[x_jx_k]=μ_jμ_k+Σ_{jk}`` for exact integration,
including covariance terms. `_cluster_θ_score` follows equation (5):

```math
\psi_\theta=\frac{W}{N}\sum_j\left[
\frac{I(A=a)}{\pi_a}\{\omega^{(1)}(Y_j-\eta_j)+u_{1j}-u_{3j}\}
+\frac{I(A=a′)}{\pi_{a′}}(u_{2j}-u_{3j})+u_{3j}\right],
```

where ``π_a=π^a(1-π)^{1-a}`` and
``ω^{(1)}=κ(M|a′,C,N)/ξ(M|a,Z,C,N)``. `_cluster_τ_score` follows equations
(6)--(7). Ratios use log densities and are evaluated only on the residual's
assigned arm. Off-arm residual terms remain exactly zero.

All members stay in their cluster's treatment-stratified fold. The ratio
estimator and cluster influence curve are

```math
\widehat\theta=\frac{\mathbb P_K\widehat\psi_\theta}{\mathbb P_KW},
\qquad
\widehat D_i=\frac{\widehat\psi_{\theta,i}-W_i\widehat\theta}{\mathbb P_KW}.
```

Linear contrasts are applied to joint influence curves before their sample
covariance is divided by ``K``. This preserves component covariance and the
individual-average normalisation. Intervals use a Student ``t_{K-1}`` critical
value as a finite-cluster convention, not an exact small-sample result.

## Optional weighted intercept fluctuation

`stabilise=true` implements Ham et al.'s Appendix I.2 adjustment. For each
validation fold and each of ``θ(0,0),θ(1,0),θ(1,1),τ``, define the weighted
outcome residual contribution ``R_i`` and residual weight ``H_i``:

```math
R_i=\frac{W_i}{N_i}\sum_j h_{ij}(Y_{ij}-\widehat\eta_{ij}),\qquad
H_i=\frac{W_i}{N_i}\sum_j h_{ij},\qquad
\widehat\epsilon=\frac{\sum_{i\in I_r}R_i}{\sum_{i\in I_r}H_i}.
```

The member weight ``h_{ij}`` includes the arm indicator, inverse randomisation
probability and mediator density ratio. Shift the outcome regression and all
compatible integrated means by the same ``ε``. Their differences are unchanged,
so the score update is directly

```math
\widehat\psi_i^{\,*}=\widehat\psi_i+(W_i-H_i)\widehat\epsilon.
```

This solves ``\sum_{i\in I_r}(R_i-H_i\widehat\epsilon)=0`` without recomputing
the integrals or changing the target. Cluster-average and individual-average
fits use different fluctuations. The adjustment changes finite-sample point
estimates; diagnostics expose `fold_fluctuations` and `fold_residual_sums`.
The ordinary trial `te` is unchanged. `stabilise=false` retains the original
reference route.

This option is separate from `ratio_cap`. Capped fits still withhold uncertainty.
Uncapped stabilised intervals require the appropriate nuisance and normalisation
rates; solving the residual equation does not establish interval validity under
misspecification or inadequate overlap. The authors omit a detailed asymptotic
justification of this adjustment in Section 6, referring to earlier work.

The corresponding R routines are
[`ap_stabilization_epsilon`, `ap_score_theta_fold` and `ap_score_tau_fold`](https://github.com/daeyoungspam-stat/CRT-interventional-mediation/blob/main/inst/engine/annals_canonical.R#L392-L475).
On 32 shared score fixtures spanning both weighting schemes, θ/τ kernels and
stabilisation settings, maximum absolute score discrepancy was
``4.09\times10^{-14}`` and fluctuation discrepancy ``3.61\times10^{-16}``.
This is score-level concordance, not validation of fitted nuisances or their
complete study. The reference file's Git blob SHA was
`f162da6239e03709ea42870c25e1f147aa9c7c3b`.

## Example

```@example cluster-mediation
using CausalMediation
include(joinpath(dirname(pathof(CausalMediation)), "..", "examples", "cluster_mediation.jl"))
cluster_example = cluster_mediation_example()
cluster_example.result.table
```

This synthetic example has a treatment-dependent confounder-mediator outcome
interaction, with different analytic ordinary and stochastic total effects.

## Inference and verification limits

Uncapped fits return `inference_status=:asymptotic`. EIF uncertainty requires
correct nuisance limits, suitable convergence rates, finite score moments and
adequate overlap; cross-fitting alone does not guarantee nominal coverage.
The general paper's multiple-robustness claims are not a guarantee for arbitrary
misspecification of this constrained parametric backend.

A finite `ratio_cap≥1` caps ratios to ``[1/cap,cap]`` for point fits and sets
`:truncated_weights`, withholding standard errors, covariance, influence curves
and limits for the whole result. Uncapped overflow or underflow throws rather
than silently changing weights.

Inspect `diagnostics.min_log_ratio`, `max_log_ratio` and `capped_fraction`.
The extrema describe active density ratios before capping. Off-arm residual
ratios are not evaluated. These diagnostics describe numerical overlap; they
do not establish the causal positivity assumption.

## Initial calibration, 13 November 2026

The CDCS driver `scripts/stress_harness/run_cluster_mediation_reference.jl`
ran 400 replications at each of 100 and 300 clusters under Julia 1.13.3, using
the example's exact Gaussian DGP and treatment-dependent ``Z\times M`` outcome
interaction. Assignment probability was 0.5, weighting was cluster-average,
cross-fitting used three folds, and ratios were uncapped. Replication ``r`` at
``K`` clusters used `MersenneTwister(10_000K + r)`.

| Effect | Coverage at K=100 | Coverage at K=300 |
|---|---:|---:|
| `ide` | 91.75% | 91.00% |
| `iie` | 95.00% | 94.25% |
| `iime` | 93.25% | 92.50% |
| `isme` | 94.50% | 92.75% |
| `stochastic_total` | 94.50% | 97.50% |
| `te` | 95.00% | 96.00% |

Coverage Monte Carlo standard errors range from 0.78 to 1.43 percentage points.
Direct-effect empirical SD was 2.824 at K=100 and 0.382 at K=300, compared with
mean reported SEs of 0.411 and 0.068. Rare extreme estimates dominate these
dispersions. For example, replication 125 at K=100 estimated `ide=-54.143`
against truth 0.56. A treated validation cluster had a control-mediator versus
treated-conditional-mediator density ratio of 5,125.656, amplifying a mean
outcome residual of -0.521.

These results establish a limitation of this **uncapped reference**, not a
comparison with the authors' stabilised implementation. The estimator is not
ready for routine inferential use. Agreement with the analytic integrals and
score equations does not establish stable fitted estimation or calibrated
intervals. Further overlap diagnostics, justified stabilisation and independent
replicated validation are needed before promoting it beyond reference use.

## Implementation comparison and ablation

The authors' canonical code uses the weighted intercept fluctuation above,
ratios bounded to ``[0.05,20]``, a selected outcome basis and common residual
covariance parameters for each fitted process across treatment arms. It derives
the parametric marginal mediator law from its fitted ``M|Z`` and ``Z`` laws.
Our Gaussian joint law already preserves that compatibility through Gaussian
marginalisation and conditioning. Pooling covariance or reducing the outcome
basis would impose additional modelling restrictions and has not been adopted
silently. Their exchangeable-density formulas offer a future optimisation;
our current dense calculations provide a more general reference.

The authors' interval convention is ``\sqrt{K/(K-3)}t_{K-3}`` applied to their
population-variance SE, whereas this reference uses sample covariance with
``t_{K-1}``. This small critical-value difference has not been used to tune away
the observed calibration problem. See their
[shared core](https://github.com/daeyoungspam-stat/CRT-interventional-mediation/blob/main/inst/engine/annals_core.R)
(`ap_fit_process`, `ap_derive_compatible_kappa`, `ap_effect_rows`;
Git blob SHA `c4c326dfc4a24b94aa4fc436527729309c985ad3`).

Using identical datasets and folds, 400 replications per cluster count gave:

| Route | IDE SD, K=100 | IDE SD, K=300 | IDE coverage, K=100 | IDE coverage, K=300 |
|---|---:|---:|---:|---:|
| Raw, uncapped | 2.824 | 0.382 | 91.75% | 91.00% |
| Intercept-adjusted, uncapped | 0.0761 | 0.0346 | 82.25% | 86.25% |
| Capped at 20 | 0.1081 | 0.0362 | Unavailable | Unavailable |
| Intercept-adjusted, capped at 20 | 0.0731 | 0.0328 | Unavailable | Unavailable |

For the uncapped adjusted route, IDE bias was 0.00137 at K=100 and -0.00161
at K=300 (Monte Carlo SEs 0.00381 and 0.00173). Point stability improved
substantially, but interval coverage remained inadequate. The adjustment is
therefore opt-in, not a claim of repaired inference. Capped routes retain the
package's unavailable-uncertainty contract; their coverage is undefined, not
zero. These results concern our own DGP and backend, not the authors' reported
simulation performance. Further work must examine the variance of the actual
stabilised estimator and validate coverage before recommending its intervals.

Centering and scaling the polynomial inputs preserved these original-scale
raw results. A regression test additionally verifies that adding ``10^6`` to
the confounder's measurement origin preserves estimates and SEs; the previous
unscaled polynomial fit failed that test.

Tests cover analytic intervention integrals, nonlinear recovery, score terms,
singleton reduction, weighting, joint covariance and invalid inputs. They do
not establish end-to-end fitted concordance with the authors' R engine.
Further calibration and that comparison remain necessary before routine use.

The existing `MediationSpec` `moc` path uses a different confounder-mediator
coupling. This API does not replace it or enable its unavailable uncertainty.
Flexible learners, the low-dimensional classification estimator, missingness
and longitudinal mediation remain outside this reference implementation.

## Hierarchical Bayesian clustered mediation

`fit_bayesian_cluster_mediation` and `bayesian_mediation_effects` provide an
experimental **parametric Bayesian** alternative for the same clustered targets
(`te`, `stochastic_total`, `ide`, `iie`, `iime`, and `isme`), reusing
`ClusterMediationSpec` (`weighting` and `π` metadata) while remaining separate
from the frequentist `run_cluster_mediation` Theorem 2 EIF reference. Rather
than forming inverse-probability and mediator density ratios, this route models
the joint law parametrically and propagates parameter uncertainty from the
posterior. It is not a semiparametric or doubly robust EIF estimator: it trades
density-ratio instability for explicit parametric mechanism forms, Gaussian
random-intercept assumptions, and declared prior specifications. Causal
identification remains `:user_assumed`.

The model factorises across three one-level random-intercept Gaussian
mechanisms:

- Post-treatment confounder ``p(Z \mid A, C, N)`` linear in ``[1, A, N, C_j, \bar{C}]``.
- Mediator ``p(M \mid A, Z, C, N)`` linear in ``[1, A, N, C_j, \bar{C}, Z_j, \bar{Z}_{-j}]``.
- Outcome ``p(Y \mid A, Z, M, C, N)`` using ``[\phi_d(x_{kj}), A_k \phi_d(x_{kj})]``,
  where ``\phi_d`` is the linear (`outcome_degree=1`) or quadratic
  (`outcome_degree=2`, default) polynomial in ``[N, C_j, \bar{C}, Z_j, \bar{Z}_{-j}, M_j, \bar{M}_{-j}]``.

Each mechanism is fitted via `CausalTargeted`'s structured variational engine
``q(\beta, b)\,q(\sigma_\varepsilon^{-2})\,q(\sigma_b^{-2})`` (Nolan, Menictas
& Wand 2020) with exact Schur-complement block algebra and exact affine design
reparameterisation. Query `cluster_mediation_design_dimensions(spec; outcome_degree)`
to construct explicit `ClusterGaussianPrior` entries for `:confounder`,
`:mediator`, and `:outcome`; no outcome-derived prior scaling is performed.

For each shared draw from the three fitted mechanism posteriors,
`bayesian_mediation_effects` integrates analytically over the outcome random
intercept and over ``(Z_k, M_k)`` using exact ``O(1)`` exchangeable-covariance
algebra (without allocating ``N_k \times N_k`` cluster matrices), returning
joint posterior means, SDs, credible intervals, MCSEs, and joint covariance
under `:cluster` or `:individual` weighting. If any mechanism fit has
`status !== :converged`, intervals and joint covariance are withheld (`NaN` /
`intervals_available = false`).

Loading `PSIS` (`using PSIS`, activating `CausalTargetedPSISExt`) enables
`correction = :psis`, which applies Pareto-smoothed importance sampling
(Vehtari et al. 2024) across the joint global parameter draws. Operational gates
(`khat < 0.7`, `ess >= 400`, and `mcse / sd <= 0.05` across non-degenerate
effects) withhold intervals and covariance (`NaN` / `intervals_available = false`)
on failure without silently falling back to raw variational intervals.

```julia
using CausalMediation, CausalTargeted, LinearAlgebra, PSIS

dims = cluster_mediation_design_dimensions(spec; outcome_degree = 2)
mkprior(p) = ClusterGaussianPrior(
    mean = zeros(p),
    covariance = 4.0 * Matrix{Float64}(I, p, p),
    residual_shape = 3.0,
    residual_scale = 2.0,
    cluster_shape = 3.0,
    cluster_scale = 1.0,
)
priors = (
    confounder = mkprior(dims.confounder),
    mediator = mkprior(dims.mediator),
    outcome = mkprior(dims.outcome),
)

fit = fit_bayesian_cluster_mediation(spec, train_df; priors, outcome_degree = 2)
res = bayesian_mediation_effects(fit, target_df; correction = :psis)
res.estimates
```

References: Ham, Cheng & Li (2026, [arXiv:2611.11903](https://arxiv.org/abs/2611.11903));
Nolan, Menictas & Wand (2020, *Streamlined computing for variational inference
with higher level random effects*, *JMLR* 21(157):1–62); Vehtari, Simpson,
Gelman, Yao & Gabry (2024, *Pareto smoothed importance sampling*, *JMLR*
25(72):1–58).
