# Length
const m = 1.0
const cm = 1e-2
const mm = 1e-3
const km = 1e3

# Time
const s = 1.0
const min = 60.0
const hr = 3600.0
const hour = 3600.0
const h = 3600.0
const day = 86400.0

# Angles
const rad = 1.0
const mrad = 1e-3
const μrad = 1e-6
const arcsec = π/648000
const arcmin = π/10800
const deg = π/180
const mdeg = 1e-3*π/180
const μdeg = 1e-6*π/180

# Weight
const kg = 1.0

# Force
const N = 1.0
const kN = 1e3

# Angular rate
const rpm = 2π/min

# Constants
const SPEED_OF_LIGHT = 299792458.0      # [m/s] Speed of light
const AU = 149597871e3                  # [m] Astronomical unit
const G0 = 9.807                        # [m/s^2] Earth gravity
const G = 6.67384e-11         	        # [m^3/kg/s^2] Universal gravitational constant
const PLANCK = 6.62607004E-34      	    # [m^2 kg/s] Planck's constant
const BOLTZMANN = 1.380649E-23          # [m^2 kg/s^2/K] Boltzmann constant
const STEFBOLTZMANN = 5.670374419e-8    # [W/m^2/K^4] Stefan-Boltzmann constant

# Conversion functions
units(u::String) = eval(Meta.parse(u))
units(val, u::String) = val.*units(u)
convertUnits(from::String, to::String) = eval(Meta.parse(from))./eval(Meta.parse(to))
convertUnits(val, from::String, to::String) = val.*convertUnits(from, to)
