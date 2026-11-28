module CausalMediationLuxExt

using CausalMediation
using Lux
using Random
using LinearAlgebra

"""
    fit_riesz_representer(X, y; hidden=16, lambda=1e-3, rng=Random.default_rng()) -> NamedTuple

Fit a lightweight single-hidden-layer `Lux` MLP Riesz-representer scaffold
(`Lux.Chain(Lux.Dense(p => hidden, tanh), Lux.Dense(hidden => 1))`).

The hidden layer weights and biases are randomly initialised via `Lux.setup`
and held fixed as a random feature map; the linear readout layer is fitted in
closed form by ridge regression with penalty `lambda` (CPU-only, no AD or
gradient training loop). GLM/SL remains the default for small-*n* book use.
"""
function fit_riesz_representer(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real};
    hidden::Int = 16,
    lambda::Real = 1e-3,
    rng::AbstractRNG = Random.default_rng(),
)
    n, p = size(X)
    length(y) == n || throw(DimensionMismatch("length(y) ($(length(y))) must match size(X, 1) ($n)"))
    hidden >= 1 || throw(ArgumentError("hidden must be >= 1; got $hidden"))
    λ = Float64(lambda)
    λ >= 0 || throw(ArgumentError("lambda must be >= 0; got $lambda"))

    hidden_layer = Lux.Dense(p => hidden, tanh)
    readout_layer = Lux.Dense(hidden => 1)
    model = Lux.Chain(hidden_layer, readout_layer)
    ps_raw, st = Lux.setup(rng, model)
    ps_h = (
        weight = Float64.(ps_raw.layer_1.weight),
        bias = Float64.(ps_raw.layer_1.bias),
    )

    Xf = Matrix{Float64}(X)
    Φ_t, st_h = Lux.apply(hidden_layer, Xf', ps_h, st.layer_1)
    Φ = Matrix{Float64}(Φ_t')
    Φ1 = hcat(ones(Float64, n), Φ)
    β = (Φ1' * Φ1 + λ * I) \ (Φ1' * Float64.(y))
    ps_fit = (
        layer_1 = ps_h,
        layer_2 = (
            weight = reshape(Float64.(β[2:end]), 1, hidden),
            bias = Float64[β[1]],
        ),
    )
    st_fit = (layer_1 = st_h, layer_2 = st.layer_2)

    return (
        kind = :lux_riesz_mlp,
        model = model,
        ps = ps_fit,
        st = st_fit,
        hidden = hidden,
        lambda = λ,
        predict = (Xnew::AbstractMatrix{<:Real} -> begin
            yhat, _ = Lux.apply(model, Matrix{Float64}(Xnew)', ps_fit, st_fit)
            return vec(Float64.(yhat))
        end),
    )
end

end # module
