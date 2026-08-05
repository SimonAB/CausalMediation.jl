# Package boundaries

**Design principles:** [DESIGN.md](DESIGN.md) · [shared](DESIGN_PRINCIPLES.md)

## CausalMediation.jl (this package)

- Mediation estimands: interventional (RI), natural, organic, controlled direct, recanting twin
- Cross-fitted EIF / one-step / TMLE / plugin nested-MC
- `moc` intermediate confounding; δ-grids; target-trial constructors
- Optional Lux Riesz-representer extension

## CausalTargeted.jl

- LMTP, sequential LMTP, Super Learner, `ShiftPolicy`, DiD, g-computation
- Soft façades that forward mediation APIs here when this package is loaded

## CausalDynamics.jl

- Graphs, `MediationQuery` / `identify`, certificates (`moc`, effect kind)

## Out of scope (for now)

- Survival / competing-risks mediation
- Default GPU deep Riesz nets
- Biological pathway registries
