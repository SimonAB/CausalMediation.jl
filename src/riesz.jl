"""Optional Riesz-representer hooks (Lux extension)."""

"""
    riesz_available() -> Bool

Return `true` when the Lux extension is loaded.
"""
function riesz_available()
    return Base.get_extension(@__MODULE__, :CausalMediationLuxExt) !== nothing
end

"""
    fit_riesz_representer(X, y; hidden=16, lambda=1e-3, rng=Random.default_rng()) -> NamedTuple

Fit a single-hidden-layer `Lux` MLP Riesz-representer scaffold (random hidden
layer via `Lux.setup` with a closed-form ridge-fitted linear readout) for
high-dimensional mediators or `moc`. Requires `using Lux` so the weakdep
extension loads. Without Lux, throws.
"""
function fit_riesz_representer(X::AbstractMatrix, y::AbstractVector; kwargs...)
    ext = Base.get_extension(@__MODULE__, :CausalMediationLuxExt)
    ext === nothing && throw(ErrorException(
        "fit_riesz_representer requires Lux.jl. Add Lux and `using Lux` to load the extension.",
    ))
    return ext.fit_riesz_representer(X, y; kwargs...)
end
