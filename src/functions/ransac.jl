# distAll = funDist(idFit) // fit using y[idFit], evaluate distance on entire y
# Robust Parameter Estimation using RANSAC
#
# The 'funDist(idFit)' function shall estimate the problem's parameters using measurements
# identified by the 'idFit' variable (e.g., y[idFit]), and shall return the distance metric
# of all measurements (inlcuding those used for the fit). The distance shall be positive, it
# is recommended to use as distance an error metric such as d = e².
#
# https://www.cse.iitb.ac.in/~ajitvr/CS763_Spring2017/RobustMethods.pdf
# https://web.ipac.caltech.edu/staff/fmasci/home/astro_refs/LeastMedianOfSquares.pdf
#
# Torr, Zisserman, MLESAC: A new robust estimator with application to estimating image geometry
# https://www.robots.ox.ac.uk/~vgg/publications/2000/Torr00/torr00.pdf
#
# Author: F. Capolupo
# European Space Agency, 2023
function ransac(
    funDist,        # Fit + distance function
    N,                      # Total number of available measurements/data points
    Nmin;                   # Minimum number of measurements required to compute a fit
    maxIter = 100,            # Maximum number of iterations
    threshold = 0.1,          # Inliers distance threshold
    verbose = true,
)

    # Init parameters and allocations
    Cbest = Inf
    dist = zeros(N, 1)
    inliers = BitArray(undef, N)
    iFit = zeros(Int, Nmin)

    # Init iterations
    for _ = 1:maxIter
        # Randomly select Nmin points to be used to fit the model
        iFit .= randperm(N)[1:Nmin]

        # Fit the model using selected Nmin random points and calculate
        # distance on entire measurement population
        dist .= abs.(funDist(iFit))     # added abs just in case (dumb user)

        # Identify inliers
        C = sum(min.(dist, threshold))

        # Update best fit if needed
        if C < Cbest
            Cbest = C
            inliers .= dist .≤ threshold
        end
    end

    if verbose
        println("RANSAC: $(sum(.!inliers)) outliers detected")
    end

    return inliers
end
