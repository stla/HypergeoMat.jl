push!(LOAD_PATH,"../src/")

using Documenter
using HypergeoMat

makedocs(
    sitename = "HypergeoMat.jl",
    authors="Stephane Laurent",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stla.github.io/HypergeoMat.jl",
        assets=String[],
    ),
    modules = [HypergeoMat],
    pages = ["Documentation"  => "index.md"],
    repo = "https://github.com/stla/HypergeoMat.jl/blob/{commit}{path}#{line}"
)

deploydocs(;
branch = "gh-pages",
    devbranch = "main",
    repo="github.com/stla/HypergeoMat.jl",
)
