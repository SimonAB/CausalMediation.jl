"""Compatibility boundary for the registered CausalTargeted 0.3.6 API."""

# CausalTargeted 0.3.7 is expected to provide the fitted schema methods used by
# the categorical-covariate path. Until that release is registered, keep the
# fallback here so estimator code does not contain version checks.
const _HAS_CAUSALTARGETED_SCHEMA = all(
    name -> isdefined(CausalTargeted, name),
    (:CovariateSchema, :fit_covariate_schema, :transform_covariates),
)

struct _CovariateSchema
    covariates::Vector{Symbol}
    backend::Any
    categorical_levels::Dict{Symbol, Vector{Any}}
end

function _fit_covariate_schema(df::AbstractDataFrame, covariates::AbstractVector{Symbol})
    requested = collect(Symbol, covariates)
    length(unique(requested)) == length(requested) || throw(ArgumentError(
        "covariates must be unique; received $(repr(requested))",
    ))
    missing_columns = [column for column in requested if !hasproperty(df, column)]
    isempty(missing_columns) || throw(ArgumentError(
        "missing requested column(s): $(join(string.(missing_columns), ", "))",
    ))
    for covariate in requested
        any(ismissing, df[!, covariate]) && throw(ArgumentError(
            "covariate :$covariate contains missing values; apply missing-data handling first",
        ))
    end

    if _HAS_CAUSALTARGETED_SCHEMA
        fit_schema = getfield(CausalTargeted, :fit_covariate_schema)
        return _CovariateSchema(requested, fit_schema(df, requested), Dict{Symbol, Vector{Any}}())
    end

    levels = Dict{Symbol, Vector{Any}}()
    for covariate in requested
        column = df[!, covariate]
        all(value -> value isa Real, column) && continue
        levels[covariate] = Any[value for value in unique(column)]
    end
    return _CovariateSchema(requested, nothing, levels)
end

function _design_matrix(
    schema::_CovariateSchema,
    df::AbstractDataFrame;
    treatment::Union{Symbol, Nothing} = nothing,
    treatment_values::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    schema.backend === nothing || return CausalTargeted.design_matrix(
        schema.backend,
        df;
        treatment = treatment,
        treatment_values = treatment_values,
    )

    n = nrow(df)
    missing_columns = [column for column in schema.covariates if !hasproperty(df, column)]
    isempty(missing_columns) || throw(ArgumentError(
        "missing requested column(s): $(join(string.(missing_columns), ", "))",
    ))
    widths = map(schema.covariates) do covariate
        haskey(schema.categorical_levels, covariate) ?
            max(length(schema.categorical_levels[covariate]) - 1, 0) : 1
    end
    X = Matrix{Float64}(undef, n, 1 + Int(treatment !== nothing) + sum(widths))
    X[:, 1] .= 1.0
    output_column = 2

    if treatment !== nothing
        hasproperty(df, treatment) || treatment_values !== nothing || throw(ArgumentError(
            "missing requested column: $treatment",
        ))
        values = treatment_values === nothing ? df[!, treatment] : treatment_values
        length(values) == n || throw(DimensionMismatch(
            "treatment values have length $(length(values)); expected $n",
        ))
        X[:, output_column] .= Float64.(values)
        output_column += 1
    end

    for (covariate, width) in zip(schema.covariates, widths)
        values = df[!, covariate]
        any(ismissing, values) && throw(ArgumentError(
            "covariate :$covariate contains missing values; apply missing-data handling first",
        ))
        if haskey(schema.categorical_levels, covariate)
            levels = schema.categorical_levels[covariate]
            unseen = unique(value for value in values if all(level -> !isequal(value, level), levels))
            isempty(unseen) || throw(ArgumentError(
                "covariate :$covariate contains unseen categorical level(s): $(join(repr.(unseen), ", "))",
            ))
            for level in Iterators.drop(levels, 1)
                X[:, output_column] .= Float64.(isequal.(values, Ref(level)))
                output_column += 1
            end
        else
            X[:, output_column] .= Float64.(values)
            output_column += width
        end
    end
    return X
end

_uses_ipcw_weights(weights::AbstractVector{<:Real}) =
    any(weight -> !isapprox(weight, 1.0; atol = 1e-12, rtol = 0.0), weights)

function _weighted_influence_summary(
    values::AbstractVector{<:Real},
    weights::AbstractVector{<:Real},
)
    length(values) == length(weights) || throw(DimensionMismatch(
        "values and weights must have the same length",
    ))
    total_weight = sum(weights)
    total_weight > 0 || throw(ArgumentError("weights must have positive total weight"))
    estimate = sum(weights .* values) / total_weight
    normalized_weights = weights ./ mean(weights)
    ic = normalized_weights .* (values .- estimate)
    return (estimate = estimate, se = std(ic) / sqrt(length(ic)), ic = Float64.(ic))
end
