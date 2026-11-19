# Bayesian Clustered Stochastic Mediation Contract

This document records the mathematical and software contract for
`fit_bayesian_cluster_mediation` and `bayesian_mediation_effects` in
`CausalMediation.jl` (Milestone 3 of `plans/cluster-effects-bayesian-inference.md`).

## 1. Relationship to `run_cluster_mediation` and `MediationSpec`

- Reuses `ClusterMediationSpec` (same columns, `weighting`, and `π` metadata).
- Remains completely separate from the frequentist `run_cluster_mediation`
  canonical-score route and does not modify `MediationSpec` or its `moc`
  sampler.
- Records causal identification assumptions (cluster randomisation, no
  cross-cluster interference, consistency, measured confounder control,
  positivity) as `:user_assumed`.

## 2. Factorised hierarchical Gaussian model

For cluster $k = 1, \dots, K$ of size $N_k$ and member $j \in \{1, \dots, N_k\}$,
let $B_{kj} = [N_k,\; C_{kj},\; \bar{C}_k]$ denote baseline features (size,
member covariates, cluster-mean covariates), and let $\bar{Z}_{k,-j}$ and
$\bar{M}_{k,-j}$ denote leave-one-out peer averages (defined as $0$ for
singleton clusters $N_k = 1$).

Three separate one-level random-intercept Gaussian mechanisms are fitted via
`CausalTargeted`'s structured variational engine
$q(\beta, b)\,q(\sigma_\varepsilon^{-2})\,q(\sigma_b^{-2})$:

1. **Post-treatment confounder** $p(Z \mid A, C, N)$:
   $$Z_{kj} = [1,\; A_k,\; B_{kj}] \beta_Z + b_{Z,k} + \varepsilon_{Z,kj}, \qquad b_{Z,k} \sim \mathcal{N}(0, \sigma_{b,Z}^2), \quad \varepsilon_{Z,kj} \sim \mathcal{N}(0, \sigma_{\varepsilon,Z}^2).$$
2. **Mediator** $p(M \mid A, Z, C, N)$:
   $$M_{kj} = [1,\; A_k,\; B_{kj},\; Z_{kj},\; \bar{Z}_{k,-j}] \beta_M + b_{M,k} + \varepsilon_{M,kj}, \qquad b_{M,k} \sim \mathcal{N}(0, \sigma_{b,M}^2), \quad \varepsilon_{M,kj} \sim \mathcal{N}(0, \sigma_{\varepsilon,M}^2).$$
3. **Outcome** $p(Y \mid A, Z, M, C, N)$:
   $$Y_{kj} = \bigl[\phi_d(x_{kj}),\; A_k \phi_d(x_{kj})\bigr] \beta_Y + b_{Y,k} + \varepsilon_{Y,kj}, \qquad x_{kj} = [B_{kj},\; Z_{kj},\; \bar{Z}_{k,-j},\; M_{kj},\; \bar{M}_{k,-j}],$$
   where $\phi_d$ is the degree-$d$ polynomial basis (`outcome_degree = 1`
   linear or `outcome_degree = 2` quadratic, matching `_cluster_polynomial`).
   The treatment-interacted block $A_k \phi_d(x_{kj})$ allows arm-specific
   interactions (including $A_k Z_{kj} M_{kj}$, so `te` and `stochastic_total`
   can differ) and supports explicit prior shrinkage on treatment interactions.

Residual innovations and cluster random intercepts are independent across the
three mechanisms conditional on their declared parents. Explicit
`ClusterGaussianPrior` objects must be supplied for `:confounder`, `:mediator`,
and `:outcome`; no outcome-dependent prior scaling is performed.

## 3. Exact Gaussian integration per posterior draw

For each joint posterior draw $\theta^{(s)} = (\beta_Z, \sigma_{\varepsilon,Z}^2, \sigma_{b,Z}^2, \beta_M, \sigma_{\varepsilon,M}^2, \sigma_{b,M}^2, \beta_Y)$
and each target cluster $k$ of size $N$:

- Confounder law under arm $a \in \{0, 1\}$:
  $\mu_Z(a) = X_Z(a)\beta_Z$, $\Sigma_Z = \sigma_{\varepsilon,Z}^2 I_N + \sigma_{b,Z}^2 J_N$.
- Conditional mediator operator $G_M$: $G_M = \gamma_{\text{own}} I_N$ for $N = 1$,
  and $G_M = \gamma_{\text{own}} I_N + \frac{\gamma_{\text{peer}}}{N - 1}(J_N - I_N)$
  for $N > 1$.
- Marginal mediator law integrating out $Z_k$:
  $$\mu_M(a) = \alpha_M(a) + G_M \mu_Z(a), \qquad \Sigma_{ZM} = \Sigma_Z G_M^\top, \qquad \Sigma_M = \Sigma_{M|Z} + G_M \Sigma_Z G_M^\top.$$
- Integrating over $b_{Y,k} \sim \mathcal{N}(0, \sigma_{b,Y}^2)$ ($\mathbb{E}[b_{Y,k}] = 0$)
  and over $(Z_k, M_k)$ analytically via exact quadratic Gaussian moments yields:
  - $\theta(a, a')$: outcome arm $a$, confounder $(\mu_Z(a), \Sigma_Z)$, independent mediator $(\mu_M(a'), \Sigma_M)$ with $\Sigma_{ZM} = 0$.
  - $\tau$: outcome arm $1$, confounder $(\mu_Z(1), \Sigma_Z)$, independent own mediator from arm $1$ and peer mediators from arm $0$ (with zero own–peer covariance).
  - $\mu_{\text{nat}}(a)$: outcome arm $a$, natural joint law $(\mu_Z(a), \Sigma_Z, \mu_M(a), \Sigma_M, \Sigma_{ZM})$.

Target averaging over clusters uses either `:cluster` weighting
($\frac{1}{K}\sum_k \frac{1}{N_k}\sum_{j=1}^{N_k} \cdot$) or `:individual`
weighting ($\frac{1}{N}\sum_k \sum_{j=1}^{N_k} \cdot$).

## 4. Draw-by-draw decomposition and failure contract

For every draw $s$:

- `ide` $= \theta(1, 0) - \theta(0, 0)$
- `iie` $= \theta(1, 1) - \theta(1, 0)$
- `stochastic_total` $= \theta(1, 1) - \theta(0, 0) = \text{ide} + \text{iie}$
- `iime` $= \tau - \theta(1, 0)$
- `isme` $= \theta(1, 1) - \tau$, with `iime + isme` $= \text{iie}$
- `te` $= \mu_{\text{nat}}(1) - \mu_{\text{nat}}(0)$

If any of the three mechanism fits fails to converge (`status !== :converged`),
`bayesian_mediation_effects` withholds credible intervals and joint covariance
(`NaN` / `intervals_available = false`).

## 5. Attribution

Causal targets follow Ham, Cheng & Li (2026), *Semiparametric causal mediation
analysis in cluster-randomized trials with post-treatment confounding*,
[arXiv:2611.11903](https://arxiv.org/abs/2611.11903). Structured Gaussian VI
follows Nolan, Menictas & Wand (2020), *JMLR* 21(18).
