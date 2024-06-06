module Units
    # Length
    m = 1.0
    cm = 1e-2
    mm = 1e-3
    km = 1e3

    # Time
    s = 1.0
    ms = 1e-3
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
    mdeg = π/180e3
    μdeg = π/180e6

    # Weight
    kg = 1.0

    # Force
    N = 1.0
    kN = 1e3
    mN = 1e-3

    # Angular rate
    rpm = 2π/min

    # Other
    ppm = 1e-6
    g = 9.80665
    g0 = g
    mg = g*1e-3
    μg = g*1e-6
    ng = g*1e-9

    # Frequency
    Hz = 1.0
    mHz = 1e-3
    kHz = 1e3

    # Constants
    SPEED_OF_LIGHT = 299792458.0        # [m/s] Speed of light
    AU = 149597871e3                    # [m] Astronomical unit
    G = 6.67384e-11         	        # [m^3/kg/s^2] Universal gravitational constant
    PLANCK = 6.62607004E-34      	    # [m^2 kg/s] Planck's constant
    BOLTZMANN = 1.380649E-23            # [m^2 kg/s^2/K] Boltzmann constant
    STEFAN_BOLTZMANN = 5.670374419e-8   # [W/m^2/K^4] Stefan-Boltzmann constant

    # Conversion functions
    export units, convertUnits
    units(u::String) = u == "%" ? 0.01 : eval(Meta.parse(u))
    convertUnits(from::String, to::String) = units(from)./units(to)

    units(val, u::String) = val.*units(u)
    convertUnits(val, from::String, to::String) = val.*convertUnits(from, to)
end
