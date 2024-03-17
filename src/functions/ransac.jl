# distAll = funDist(idFit) // fit using y[idFit], evaluate distance on entire y
# Robust Parameter Estimation using RANSAC
#
# https://www.cse.iitb.ac.in/~ajitvr/CS763_Spring2017/RobustMethods.pdf
# https://web.ipac.caltech.edu/staff/fmasci/home/astro_refs/LeastMedianOfSquares.pdf
#
# Author: F. Capolupo
# European Space Agency, 2023
function ransac(funDist, N, Nmin; maxIter=100, threshold=0.1)

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
        inl = dist .≤ threshold
        nIn = sum(inl)
        if nIn > nBest
            nBest = nIn
            inliers .= copy(inl)
        end
    end

    return inliers
end
