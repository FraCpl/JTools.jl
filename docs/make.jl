using JTools
using Documenter

DocMeta.setdocmeta!(JTools, :DocTestSetup, :(using JTools); recursive=true)

makedocs(;
    modules=[JTools],
    authors="F. Capolupo",
    repo="https://github.com/FraCpl/JTools.jl/blob/{commit}{path}#{line}",
    sitename="JTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", canonical="https://FraCpl.github.io/JTools.jl", edit_link="master", assets=String[]
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/FraCpl/JTools.jl", devbranch="master")
