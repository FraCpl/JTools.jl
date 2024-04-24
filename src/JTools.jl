module JTools

#using Plots
using ForwardDiff
using FiniteDiff
using GLMakie
using FileIO
using GeometryBasics
using LinearAlgebra
using Random
using SpecialFunctions
using FFTW

export ecdf, goldenSectionSearch, rootFinder, ode38, trapz, cumtrapz
export dm2dv, dm2dt, dv2dm, dv2dt
export theme_fra, plotMoon!, plotEarth!, axisoff!, plotframe!, transformModel, sensorFovModel, cuboidModel
export plotEcdf, plotEcdf!, plotBox, plotBox!, lines3!
export nicholsgrid
export crossmat, logrange, mag2db, db2mag, unwrap!
export lsq
export ransac
export montecarloNsim, montecarloConfidence
export gradientDescent
export psd, psd2var
include("functions/trapz.jl")
include("functions/ecdf.jl")
include("functions/goldenSectionSearch.jl")
include("functions/rootFinder.jl")
include("functions/ode38.jl")
include("functions/rocketEqs.jl")
include("functions/makieTools.jl")
include("functions/math.jl")
include("functions/lsq.jl")
include("functions/ransac.jl")
include("functions/montecarlo.jl")
include("functions/gradientDescent.jl")
include("functions/psd.jl")

include("functions/units.jl")
using .Units
export units, convertUnits

end
