Base.@kwdef struct UNITSDATA
    # Length
    m = 1.0
    cm = 1e-2
    mm = 1e-3
    km = 1e3

    # Time
    s = 1.0
    min = 60.0
    hr = 3600.0
    hour = 3600.0
    h = 3600.0
    day = 86400.0

    # Angles
    rad = 1.0
    mrad = 1e-3
    μrad = 1e-6
    arcsec = π/648000
    arcmin = π/10800
    deg = π/180
    mdeg = 1e-3*π/180
    μdeg = 1e-6*π/180

    # Weight
    kg = 1.0

    # Force
    N = 1.0
    kN = 1e3

    # Angular rate
    rpm = 2π/min

    # Constants
    SPEED_OF_LIGHT = 299792458.0        # [m/s] Speed of light
    AU = 149597871e3                    # [m] Astronomical unit
    G0 = 9.807                          # [m/s^2] Earth gravity
    G = 6.67384e-11         	        # [m^3/kg/s^2] Universal gravitational constant
    PLANCK = 6.62607004E-34      	    # [m^2 kg/s] Planck's constant
    BOLTZMANN = 1.380649E-23            # [m^2 kg/s^2/K] Boltzmann constant
    STEFBOLTZMANN = 5.670374419e-8      # [W/m^2/K^4] Stefan-Boltzmann constant
end

_UNITSDATA_ = UNITSDATA()

# Conversion functions
units(u::String) = eval(Meta.parse("_UNITSDATA_."*u))
units(val, u::String) = val.*units(u)
convertUnits(from::String, to::String) = eval(Meta.parse("_UNITSDATA_."*from))./eval(Meta.parse("_UNITSDATA_."*to))
convertUnits(val, from::String, to::String) = val.*convertUnits(from, to)
