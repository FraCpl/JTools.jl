# This function computes the minimum number of simulations needed to verify
# a probability of success P (e.g., 0.9973) with confidence level L (e.g.
# 0.95), knowing the number of failures Nfail.
#
# References:
# [1] ECSS-E-ST-60-20C Rev.1, Section E-1
function montecarloNsim(P::Float64, L::Float64, Nfail::Int)
    for Nsim in Nfail:10000
        if mcGetConfidence(P, Nsim, Nfail) ≥ L
            return Nsim
        end
    end
end

# Confidence level of the probability of success P given Nsim simulations and Nfail failures
montecarloConfidence(P::Float64, Nsim::Int, Nfail::Int) = beta_inc(Nfail + 1, Nsim - Nfail + 1, 1.0 - P)[1]
