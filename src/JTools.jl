module JTools

#using Plots
using ForwardDiff
using FiniteDiff
using GLMakie
using FileIO
using GeometryBasics
using LinearAlgebra
using Random

export ecdf, goldenSectionSearch, rootFinder, ode38, trapz, cumtrapz
export dm2dv, dm2dt, dv2dm, dv2dt, units, convertUnits
export theme_fra, plotMoon!, axisoff!, plotframe!, transformModel, sensorFovModel, cuboidModel
export nicholsgrid
export crossmat, logrange, mag2db, db2mag, unwrap!
export lsq
export ransac
include("functions/trapz.jl")
include("functions/ecdf.jl")
include("functions/goldenSectionSearch.jl")
include("functions/rootFinder.jl")
include("functions/ode38.jl")
include("functions/rocketEqs.jl")
include("functions/units.jl")
include("functions/makieTools.jl")
include("functions/math.jl")
include("functions/lsq.jl")
include("functions/ransac.jl")

end
