# References

Bibliographic keys match the CDCS book file `references.bib` where possible, so
chapters and package docs stay aligned. Prefer DOIs when citing externally.

## Mediation (natural, interventional, stochastic)

- Robins, J. M., & Greenland, S. (1992). Identifiability and exchangeability for direct and indirect effects. *Epidemiology*, *3*(2), 143–155. — key `robins1992estimation`

- Pearl, J. (2001). Direct and indirect effects. In *UAI*. — key `pearl2001direct`

- VanderWeele, T. J. (2015). *Explanation in Causal Inference: Methods for Mediation and Interaction*. Oxford University Press. — key `vanderweele2015explanation`

- Vansteelandt, S., & Daniel, R. M. (2017). Interventional effects for mediation analysis with multiple mediators. *Epidemiology*, *28*(2), 258–265. [doi:10.1097/EDE.0000000000000596](https://doi.org/10.1097/EDE.0000000000000596) — key `vansteelandt2017interventional`

- Díaz, I., & Hejazi, N. S. (2020). Causal mediation analysis for stochastic interventions. *Journal of the Royal Statistical Society: Series B*, *82*(3), 661–683. [doi:10.1111/rssb.12362](https://doi.org/10.1111/rssb.12362) — key `diaz2020mediation`

- Hejazi, N. S., Rudolph, K. E., van der Laan, M. J., & Díaz, I. (2023). Nonparametric causal mediation analysis for stochastic interventional (in)direct effects. *Biostatistics*, *24*(3), 686–707. [doi:10.1093/biostatistics/kxac002](https://doi.org/10.1093/biostatistics/kxac002) — key `hejazi2023stochastic`

- Liu, R., Williams, N. T., Rudolph, K. E., & Díaz, I. (2024). General targeted machine learning for modern causal mediation analysis. arXiv:2408.14620. [doi:10.48550/arXiv.2408.14620](https://doi.org/10.48550/arXiv.2408.14620) — key `liu2024mediation`

- Liu, R., Williams, N. T., Rudolph, K. E., & Díaz, I. (2025). crumble: A comprehensive framework for modern causal mediation analysis with intermediate confounding. arXiv:2604.09902. [doi:10.48550/arXiv.2604.09902](https://doi.org/10.48550/arXiv.2604.09902) — key `liu2025crumble`

## Organic effects and recanting twins

- Lok, J. J. (2015). Organic direct and indirect effects with multiple mediators. *Statistics in Medicine* (and related Lok papers). — cite via book bibliography when key is present

- Vo, T. T., & Díaz, I. (and related Vo–Díaz work on recanting twins / path-specific effects). — see CDCS `references.bib` for the keyed entry used in the book

## Modified treatment policies (shared with CausalTargeted)

- Díaz Muñoz, I., & van der Laan, M. J. (2012). Population intervention causal effects based on stochastic interventions. *Biometrics*, *68*(2), 541–549. [doi:10.1111/j.1541-0420.2011.01685.x](https://doi.org/10.1111/j.1541-0420.2011.01685.x) — key `diaz2012stochastic`

- Díaz, I., Williams, N., Hoffman, K. L., & Schenck, E. J. (2023). Nonparametric causal effects based on longitudinal modified treatment policies. *Journal of the American Statistical Association*, *118*(542), 846–857. [doi:10.1080/01621459.2021.1955691](https://doi.org/10.1080/01621459.2021.1955691) — key `diaz2023lmtp`

## Targeted learning and Super Learner

- van der Laan, M. J., & Rubin, D. (2006). Targeted maximum likelihood learning. *The International Journal of Biostatistics*, *2*(1). — key `vanderlaan2006targeted`

- van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). Super learner. *Statistical Applications in Genetics and Molecular Biology*, *6*(1). — key `vanderlaan2007super`

- van der Laan, M. J., & Rose, S. (2011). *Targeted Learning*. Springer. — key `vanderlaan2011targeted`

- Zheng, W., & van der Laan, M. J. (2011). Cross-validated targeted minimum-loss-based estimation. In van der Laan & Rose (2011). — key `zheng2011crossfitting`

## Identification and target trials (upstream)

- Pearl, J. (2009). *Causality* (2nd ed.). Cambridge University Press. — key `pearl2009causality`

- Hernán, M. A., & Robins, J. M. (2020). *Causal Inference: What If*. Chapman & Hall/CRC. — key `hernan2020causal`

## Related software

- R packages [`crumble`](https://cran.r-project.org/package=crumble), [`medoutcon`](https://cran.r-project.org/package=medoutcon), [`medRCT`](https://cran.r-project.org/package=medRCT), [`lmtp`](https://cran.r-project.org/package=lmtp)

- Python: [Ananke](https://github.com/UH-CAnD3/ananke)

- Julia stack: [CausalDynamics.jl](https://github.com/SimonAB/CausalDynamics.jl), [CausalTargeted.jl](https://github.com/SimonAB/CausalTargeted.jl), [DAGMakie.jl](https://github.com/SimonAB/DAGMakie.jl)
