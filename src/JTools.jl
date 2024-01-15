module JTools

using Plots
using ForwardDiff
using FiniteDiff

export ecdf, goldenSectionSearch, rootFinder, ode38, trapz, cumtrapz
export dm2dv, dm2dt, dv2dm, dv2dt
include("functions/trapz.jl")
include("functions/ecdf.jl")
include("functions/goldenSectionSearch.jl")
include("functions/rootFinder.jl")
include("functions/ode38.jl")
include("functions/rocketEqs.jl")

end
