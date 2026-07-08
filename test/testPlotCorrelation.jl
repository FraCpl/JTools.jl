using LinearAlgebra
using JTools


data = (A=randn(1000), B=randn(1000), C=randn(1000), D=rand(1000))
isValid = data.A .> 0.5
clr = [iv ? :green : :red for iv in isValid]
plotCorrelation(data; color=clr)
