# distAll = funDist(idFit) // fit using y[idFit], evaluate distance on entire y
# Robust Parameter Estimation using RANSAC
#
# The 'funDist(idFit)' function shall estimate the problem's parameters using measurements
# identified by the 'idFit' variable (e.g., y[idFit]), and shall return the distance metric
# of all measurements (inlcuding those used for the fit).
#
# https://www.cse.iitb.ac.in/~ajitvr/CS763_Spring2017/RobustMethods.pdf
# https://web.ipac.caltech.edu/staff/fmasci/home/astro_refs/LeastMedianOfSquares.pdf
#
# Author: F. Capolupo
# European Space Agency, 2023
function ransac(funDist,        # Fit + distance function
        N,                      # Total number of available measurements/data points
        Nmin;                   # Minimum number of measurements required to compute a fit
        maxIter=100,            # Maximum number of iterations
        threshold=0.1,          # Inliers distance threshold
    )

    # Init parameters and allocations
    nBest = 0
    dist = zeros(N, 1)
    inliers = BitArray(undef, N)

    # Init iterations
    for _ in 1:maxIter
        # Randomly select Nmin points to be used to fit the model
        iFit = randperm(N)[1:Nmin]

        # Fit the model using selected Nmin random points and calculate
        # distance on entire measurement population
        dist .= abs.(funDist(iFit))     # added abs just in case (dumb user)

        # Identify inliers
        inl = dist .≤ threshold         # Inliers for current fit
        nIn = sum(inl)                  # Number of inliers
        if nIn > nBest
            # Update best fit if needed
            nBest = nIn
            inliers .= copy(inl)
        end
    end

    println("RANSAC: $(N - nBest) outliers detected")

    return inliers
end
