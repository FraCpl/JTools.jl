"""
    Nsim = montecarloNsim(P, L, Nfail)

Compute the minimum number of simulations `Nsim` needed to verify a probability of success
`P` with confidence level `L`, knowing the number of failed simulations `Nfail`. If the
number of simulations exceeds 100_000, the function returns -1.

### Examples
```julia-repl
julia> montecarloNsim(0.997, 0.95, 0)
997
```

### References
[1] ECSS-E-ST-60-20C Rev.1, Section E-1
"""
function montecarloNsim(P::Float64, L::Float64, Nfail::Int = 0)
    for Nsim = Nfail:100_000
        if montecarloConfidence(P, Nsim, Nfail) ≥ L
            return Nsim
        end
    end
    return -1
end

"""
    L = montecarloConfidence(P, Nsim, Nfail)

Compute the confidence level `L` of the probability of success `P` given `Nfail` failed
simulations over `Nsim` simulations.

### Examples
```julia-repl
julia> montecarloConfidence(0.997, 1000, 0)
0.950585606425314
```

### References
[1] ECSS‐E‐ST‐60‐20C Rev. 1 page 77
"""
montecarloConfidence(P::Float64, Nsim::Int, Nfail::Int) =
    beta_inc(Nfail + 1, Nsim - Nfail + 1, 1.0 - P)[1]
