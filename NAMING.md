# Naming

Canonical Julia names (no R brand in the public API):

| Julia | R `crumble` / literature |
|-------|---------------------------|
| `:interventional` / `InterventionalMediation` | `effect="RI"` |
| `:natural` / `NaturalMediation` | `effect="N"` |
| `:organic` / `OrganicMediation` | `effect="O"` |
| `:recanting_twin` / `RecantingTwinMediation` | `effect="RT"` |
| `moc` | intermediate confounders (`Z`) |
| `DiscreteTreatmentPolicy` on both arms | factor-`A` recode (interventional only) |
| `run_mediation` / `run_mediation_grid` | `crumble(...)` |

CausalTargeted no longer exports legacy `run_crumble_*` mediation façades; CausalMediation.jl is the primary public mediation surface.
