# Need this because GTOC 11 uses a different definition of an AU than the standard
@unit AU "AU" AstronomicalUnit 149_597_870_691*m false

# Default units
const DEFAULT_TIME_UNIT = yr
const DEFAULT_DISTANCE_UNIT = AU
const DEFAULT_MASS_UNIT = kg
const DEFAULT_ANGLE_UNIT = °

# Solver constnats
const DEFAULT_ALG = ARKODE(Explicit(), etable=FEHLBERG_13_7_8)
const DEFAULT_SIM_ARGS = (adaptive=false, dt=Float64(ustrip(DEFAULT_TIME_UNIT(1d))))

# Problem constants
const Γ = 1e-4(m/s^2) |> DEFAULT_DISTANCE_UNIT/DEFAULT_TIME_UNIT^2                    # Constant thrust acceleration
const μ = 1.32712440018e11(km^3/s^2) |> DEFAULT_DISTANCE_UNIT^3/DEFAULT_TIME_UNIT^2   # Sun gravitational parameter
const α = 6e-9(1/s) |> DEFAULT_TIME_UNIT^-1                                           # Mass proportion coefficient
