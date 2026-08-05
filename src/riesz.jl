"""Optional Riesz-representer hooks (Lux extension)."""

"""
    riesz_available() -> Bool

Return `true` when the Lux extension is loaded.
"""
function riesz_available()
    return Base.get_extension(@__MODULE__, :CausalMediationLuxExt) !== nothing
end

"""
    fit_riesz_representer(X, y; kwargs...) -> Any

Fit a Riesz representer for high-dimensional mediators / moc.
Requires `using Lux` so the weakdep extension loads. Without Lux, throws.
"""
function fit_riesz_representer(X::AbstractMatrix, y::AbstractVector; kwargs...)
    ext = Base.get_extension(@__MODULE__, :CausalMediationLuxExt)
    ext === nothing && throw(ErrorException(
        "fit_riesz_representer requires Lux.jl. Add Lux and `using Lux` to load the extension.",
    ))
    return ext.fit_riesz_representer(X, y; kwargs...)
end
