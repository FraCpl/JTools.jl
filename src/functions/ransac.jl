# RANSAC (MSAC) algorithm
#
# distAll = funDist(idFit) // fit using y[idFit], evaluate distance on entire y
# Robust Parameter Estimation using RANSAC
#
# The 'funDist(idFit)' function shall estimate the problem's parameters using data
# identified by the 'idFit' variable (e.g., y[idFit]), and shall return the distance metric
# of all measurements (inlcuding those not used for the fit). The distance shall be *positive*,
# it is recommended to use as distance an error metric such as d = e².
#
# References:
# [1] https://www.cse.iitb.ac.in/~ajitvr/CS763_Spring2017/RobustMethods.pdf
# [2] https://web.ipac.caltech.edu/staff/fmasci/home/astro_refs/LeastMedianOfSquares.pdf
# [3] Torr, Zisserman, MLESAC: A new robust estimator with application to estimating image
#     geometry, https://www.robots.ox.ac.uk/~vgg/publications/2000/Torr00/torr00.pdf
#
# Author: F. Capolupo
# European Space Agency, 2023
function ransac(
    funDist,                # Fit + distance function
    N::Int,                 # Total number of available measurements/data points
    Nmin::Int;              # Minimum number of measurements required to compute a fit
    maxIter::Int=100,       # Maximum number of iterations
    threshold=0.1,          # Inliers distance threshold
    verbose::Bool=true,
)

    # Init parameters and allocations
    Cbest = Inf
    dist = zeros(N)
    inliers = BitVector(undef, N)
    iFit = zeros(Int, Nmin)
    nPerm = zeros(Int, N)

    # Init iterations
    for _ in 1:maxIter
        # Randomly select Nmin points to be used to fit the model
        randperm!(nPerm)
        @inbounds for i in eachindex(iFit)
            iFit[i] = nPerm[i]
        end

        # Fit the model using selected Nmin random points and calculate
        # distance on entire measurement population
        dist .= funDist(iFit)

        # Identify inliers
        C = 0.0
        @inbounds for i in eachindex(dist)
            dist[i] = abs(dist[i])          # added abs just in case (dumb user)
            C += min(dist[i], threshold)    # MSAC (Slide 42 of [1])
        end

        # Update best fit if needed
        if C < Cbest
            Cbest = C
            @inbounds for i in eachindex(dist)
                inliers[i] = dist[i] ≤ threshold
            end
        end
    end

    verbose && println("RANSAC: $(sum(.!inliers)) outliers detected")

    return inliers
end
