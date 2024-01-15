# Mass-ΔV conversion functions
dv2dm(Δv, m0, Isp) = m0*(1 - exp(-Δv/9.807/Isp))
dv2dt(Δv, m0, Isp, T) = dv2dm(Δv, m0, Isp)*Isp*9.807/T
dm2dv(Δm, m0, Isp) = Isp*9.807*log(m0/(m0 - Δm))
dm2dt(Δm, m0, Isp, T) = dv2dt(dm2dv(Δm, m0, Isp), m0, Isp, T)
