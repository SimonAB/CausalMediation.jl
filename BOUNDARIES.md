# Package boundaries

**Design principles:** [DESIGN.md](DESIGN.md) · [shared](DESIGN_PRINCIPLES.md)

## CausalMediation.jl (this package)

- Mediation estimands: interventional (RI), natural, organic, controlled direct, recanting twin
- Cross-fitted EIF / one-step / TMLE / plugin nested-MC
- Numeric MTP (`ShiftPolicy`) and factor recodes (`DiscreteTreatmentPolicy`) for interventional TE/NDE/NIE (continuous `M`, empty `moc` on the factor path)
- `moc` intermediate confounding; δ-grids; target-trial constructors
- Optional Lux Riesz-representer extension

## CausalTargeted.jl

- LMTP, sequential LMTP, Super Learner, `ShiftPolicy`, DiD, g-computation
- `DiscreteTreatmentPolicy` / sequential factor recodes (LMTP; mediation EIF stays in CausalMediation)
- Soft façades that forward mediation APIs here when this package is loaded

## CausalDynamics.jl

- Graphs, `MediationQuery` / `identify`, certificates (`moc`, effect kind)

## Out of scope (for now)

- Raw high-dimensional tensors as mediators (spectra, images): compress first via
  CausalDynamics `RepresentationSpec` / `encode_to_panel`, then pass **code**
  columns to `MediationSpec`; see [Methods](docs/src/methods.md#high-dimensional-mediators-via-codes)
- Categorical mediators (`g(M|A,W)` multinomial)
- Intermediate confounding (`moc`) on the factor-`A` path
- Survival / competing-risks mediation
- Sequential / longitudinal mediation (`A_t`, `M_t`)
- Default GPU deep Riesz nets
- Biological pathway registries
- Deep generative mechanisms with non-additive (encoder) abduction on raw
  image/tensor nodes; prefer Phase 1 codes + CausalDynamics Phase 2b
  `:generative` additive L3 on those codes
