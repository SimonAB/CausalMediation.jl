# CausalMediation.jl — design principles

This package is the **mediation estimation layer**: interventional, natural, organic,
controlled direct, and recanting-twin effects with cross-fitted EIF / one-step / TMLE.

**Shared principles:** [DESIGN_PRINCIPLES.md](DESIGN_PRINCIPLES.md)  
**Boundaries:** [BOUNDARIES.md](BOUNDARIES.md)

## Role in the stack

```
IdentificationResult  →  plan_mediation / run_mediation  →  TE/NDE/NIE (+ path terms)
         ↑                        ↑
   CausalDynamics          CausalTargeted (SL, ShiftPolicy, folds)
```

## Package-specific principles

- Do **not** duplicate Super Learner; call CausalTargeted.
- Public effect names are Julia symbols (`:interventional`, `:natural`, …), not R brand names.
- Intermediate confounders are first-class (`moc` on `MediationSpec` and `MediationQuery`).
- Natural effects refuse nonempty `moc` (shared gate with `identify`).
