using Documenter
using CausalMediation

makedocs(
    sitename = "CausalMediation.jl",
    authors = "Simon A. Babayan",
    modules = [CausalMediation],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://simonab.github.io/CausalMediation.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "Naming" => "naming.md",
    ],
    checkdocs = :none,
    warnonly = [:missing_docs, :cross_references],
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/SimonAB/CausalMediation.jl.git",
        devbranch = "main",
        push_preview = true,
    )
end
