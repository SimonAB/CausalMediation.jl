# Naming

```@eval
# Include package NAMING.md body via Documenter pages; keep in sync.
```

Canonical Julia names (no R brand in the public API):

| Julia | R `crumble` / literature |
|-------|---------------------------|
| `:interventional` / `InterventionalMediation` | `effect="RI"` |
| `:natural` / `NaturalMediation` | `effect="N"` |
| `:organic` / `OrganicMediation` | `effect="O"` |
| `:recanting_twin` / `RecantingTwinMediation` | `effect="RT"` |
| `moc` | intermediate confounders (`Z`) |
| `run_mediation` / `run_mediation_grid` | `crumble(...)` |

Legacy CausalTargeted aliases `run_crumble_*` remain soft-deprecated façades.
