# Author: F. Capolupo
# European Space Agency, 2024
module JTools

using ForwardDiff
using FiniteDiff
using GLMakie
using FileIO
using GeometryBasics
using LinearAlgebra
using Random
using SpecialFunctions
using FFTW
using Printf
using Statistics
using Convex, ECOS          # RCS Conf
# using Delaunay              # RCS Conf Analysis [REMOVED --> This needs PyCall, not easy to install]

export trapz, cumtrapz
include("functions/trapz.jl")

export ecdf, interpfEcdf, interpxEcdf
include("functions/ecdf.jl")

export goldenSectionSearch
include("functions/goldenSectionSearch.jl")

export rootFinder
include("functions/rootFinder.jl")

export rootFinder!
include("functions/rootFinderInPlace.jl")

export bisection
include("functions/bisection.jl")

export ode38, ode38!
include("functions/ode38.jl")

export dm2dv, dm2dt, dv2dm, dv2dt
include("functions/rocketEqs.jl")

export theme_fra, plotMoon!, plotEarth!, axisoff!, plotframe!, transformModel, sensorFovModel
export cuboidModel, plotCuboid!, plotCone!, plotCylinder!, plotEcdf, plotEcdf!
export plotBox, plotBox!, lines3!, scatter3!, multilines!, plotSphere!, plotSensorFov!, mergeMesh
export plotCorrelation, plotCorrelation!, plotEllipse!
include("functions/makieTools.jl")

export crossmat, unwrap!, modd, isMultiple, signum #, logrange
export polyfit, polyval, interp1, mul3x1!, mul3x3!, evalpoly!
include("functions/math.jl")

export lsq, lsqWeight
include("functions/lsq.jl")

export ransac
include("functions/ransac.jl")

export montecarloNsim, montecarloConfidence
include("functions/montecarlo.jl")

export gradientDescent, Grad, RMSProp, Momentum, Adam
include("functions/gradientDescent.jl")

export psd, psd2var
include("functions/psd.jl")

export rcsMixMatrix, rcsAllocation, rcsAnalysis, rcsEnvelope, plotRcs, plotRcs!
export rcsAllocationSimplex, rcsAllocationSimplex!, RcsAllocator, analyzeRcsConf, analyzeRcsConf2D
include("functions/rcsTools.jl")

export mesh2obj, ObjModel, readObj, writeObj, grid2mesh, grid2obj
include("functions/obj.jl")

export includedir, readdirext
include("functions/utils.jl")

include("functions/units.jl")
using .Units
export units, convertUnits

end
