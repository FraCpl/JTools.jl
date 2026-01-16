using JTools

f1(x) = (x - 2.5)^2
xmin = goldenSectionSearch(f1, -10.0, 10.0)
@show xmin - 2.5 # Should be 2.5
