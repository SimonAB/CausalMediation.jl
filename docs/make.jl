using Documenter
using CausalMediation
using CausalTargeted
using CausalDynamics
using Graphs
using CairoMakie

# Prefer PNG MIME so Documenter writes figure files instead of huge inline HTML
# (same convention as CausalTargeted.jl / CausalDynamics.jl docs).
CairoMakie.activate!(type = "png")
CairoMakie.enable_only_mime!("png")

makedocs(
    sitename = "CausalMediation.jl",
    authors = "Simon A. Babayan",
    modules = [CausalMediation],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://simonab.github.io/CausalMediation.jl",
        assets = String[],
        example_size_threshold = 0,
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "Comparison" => "comparison.md",
        "Methods and literature" => "methods.md",
        "Naming" => "naming.md",
        "API overview" => "api.md",
        "References" => "references.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/SimonAB/CausalMediation.jl.git",
        devbranch = "main",
        push_preview = true,
    )
end
