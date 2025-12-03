# Verified with DSP.periodogram(w; fs=1/gyro.Ts, onesided=true)
function psd(x, Ts; twosided=false)
    # % http://www.issibern.ch/forads/sr-001-01.pdf
    # https://www.mathworks.com/help/signal/ug/power-spectral-density-estimates-using-fft.html
    #
    N = length(x)
    N -= N%2 != 0                   # Make it even
    n = round(Int, N/2)
    Y = abs.(fft(x[1:N]))
    f = 1/Ts/N*(0:n)
    Gx = 2Ts/N*(Y[1:(n + 1)]) .^ 2        # One-sided PSD
    if twosided
        Gx ./= 2.0
    end
    return f, Gx
end

# f only positive
function psd2var(f, Gx; twosided=false)
    varx = cumtrapz(f, Gx)
    return twosided ? 2varx : varx
end
