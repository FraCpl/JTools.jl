# Compute Empirical Cumulative Distribution Function F(y) of input data 'y'.
#
# Author: F. Capolupo
# European Space Agency, 2021
#
# TEST:  ecdf(abs.(randn(100000)) + 2.0*rand(100000); doPlot=true)
function ecdf(y::Vector; doPlot::Bool=false)

    # Check inputs
    x = copy(y)
    filter!(x -> ~(isnan.(x) .| isinf.(x)), x)   # Remove NaN and Inf entries

    if isempty(x)
        return NaN, NaN
    end

    # Compute Empirical F(x)
    sort!(x)
    n = length(x)
    F = (1:n)./n

    # Plot
    if doPlot
        set_theme!(theme_fra(true))
        f = Figure(resolution=(750, 500))
        ax = Axis(f[1, 1]; limits=(minimum(x), maximum(x), 0.0, 1.05), xlabel="x", ylabel="P(x), F(x)")

        hist!(ax, x; bins=100, normalization=:pdf, lab="", alpha=0.3, strokewidth=0.0, color=:gray)
        stairs!(x, F; color=:white)
        display(f)
        #=
        display(
            plot!(
                x, F;
                linewidth=2,
                color=:white,
                linetype=:steppre,
                ylims=(0.0, 1.05),
                xlabel="x",
                ylabel="P(x), F(x)",
                lab="",
                ticks=:native,
                size=(750, 500),
                framestyle=:box,
                xlim=(minimum(x), maximum(x)),
                bg=RGB(40/255, 44/255, 52/255),
                fg=RGB(0.7,0.7,0.7),
                margin =0.3Plots.cm
            )
        )
        =#
        return x, F, f
    end

    return x, F
end
