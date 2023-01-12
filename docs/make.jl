using PoissonSolvers
using Documenter

ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(PoissonSolvers, :DocTestSetup, :(using PoissonSolvers); recursive=true)

makedocs(;
    modules=[PoissonSolvers],
    authors="Michael Kraus",
    repo="https://github.com/JuliaPlasma/PoissonSolvers.jl/blob/{commit}{path}#{line}",
    sitename="PoissonSolvers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaPlasma.github.io/PoissonSolvers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => "poisson.md",
        "Library" => "library.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaPlasma/PoissonSolvers.jl",
    devbranch="main",
)
