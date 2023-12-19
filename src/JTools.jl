module JTools

using Plots
using ForwardDiff

export ecdf, goldenSectionSearch, rootFinder
include("functions/ecdf.jl")
include("functions/goldenSectionSearch.jl")
include("functions/rootFinder.jl")

end
