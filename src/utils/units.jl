#Dist
const _convSIdist = Dict("m" => 1.0, 
                        "km"=>1000.0, 
                        "nmi" => nmi_to_m,
                        "ft"=> ft_to_m, 
                        "in"=> in_to_m)

function convertDist(value::Float64, units_in::AbstractString="m", units_out="m")
    dict = _convSIdist
    return value * dict[units_in]/dict[units_out]
end

"""
    parse_unit(str::AbstractString)

Uses a regex match to split a string with magnitude and units and return them.

# Examples
```julia-repl
julia> parse_unit_string("3.2 ft")
(3.2, "ft")

julia> parse_unit_string("8m2")
(8.0, "m2")

julia> parse_unit_string("5e3  kN")
(5000.0, "kN")

julia> parse_unit_string("-32e2 W")
(-3200.0, "W")
```
"""
function parse_unit(str::AbstractString)
    m = match(r"(\+?\-?\d+\.?e?\d*)\s*(\w*\/?\w*)", str)
    mag = parse(Float64, m[1])
    units = m[2]
    if units == ""
        return mag
    else
        return mag, units
    end
end

#handle floats/ints directly, assuming SI input
function parse_unit(input::Real)
    return float.(input)
end

function get_unit_dim(str::AbstractString)
    m = match(r"([^\d]+)(\d*)", str)
    unit = m[1]
    if m[2] == ""
        dim = 1.0
    else
        dim = parse(Float64, m[2])
    end
    return unit, dim
end

function convertArea(value::Float64, units_in::AbstractString="m2", units_out="m2")
    dict = _convSIdist
    baseunit_in, exp_in = get_unit_dim(units_in)
    baseunit_out, exp_out = get_unit_dim(units_out)

    if exp_out != exp_in
        error("Units don't have the same dimensions")
    elseif exp_in != 2.0
        error("Area needs to be units of length squared (exponent = 2)")
    end

    return value * (dict[baseunit_in]^exp_in)/(dict[baseunit_out]^exp_out)
end

function convertVolume(value::Float64, units_in::AbstractString="m3", units_out="m3")
    dict = _convSIdist
    baseunit_in, exp_in = get_unit_dim(units_in)
    baseunit_out, exp_out = get_unit_dim(units_out)

    if exp_out != exp_in
        error("Units don't have the same dimensions")
    elseif exp_in != 3.0
        error("Volume needs to be units of length cubed (exponent = 3)")
    end

    return value * (dict[baseunit_in]^exp_in)/(dict[baseunit_out]^exp_out)
end

#Mass
const _convSImass = Dict("kg"=>1.0,
                         "g"=>1e-3,
                         "lbm"=>1.0/2.20462)

function convertMass(value::Float64, units_in::AbstractString="kg", units_out="kg")
    dict = _convSImass
    return value * dict[units_in]/dict[units_out]
end
function convertDensity(value::Float64, units_in:: AbstractString="kg/m3", units_out="kg/m3")
    massdict = _convSImass
    lendict = _convSIdist
    massunitin, lenunitin = split(units_in, "/")
    massunitout, lenunitout = split(units_out, "/")

    return value * (massdict[massunitin]/massdict[massunitout]) /convertVolume(1.0, lenunitin, lenunitout)
end
#Force
const _convSIforce = Dict("N" => 1.0,
                          "kN" => 1000.0,
                          "lbf" => lbf_to_N)       

function convertForce(value::Float64, units_in::AbstractString="N", units_out="N")
    dict = _convSIforce
    return value * dict[units_in]/dict[units_out]
end

const _convSIpressure = Dict("Pa" => 1.0,
                            "atm" => 101325.0,
                            "lbf/in2" => lbf_to_N/in_to_m^2,
                            "psi" => lbf_to_N/in_to_m^2,
                            "lbf/ft2" => lbf_to_N/ft_to_m^2)

function convertPressure(value::Float64, units_in::AbstractString="Pa", units_out="Pa")
    dict = _convSIpressure
    return value * dict[units_in]/dict[units_out]
end
#Speed
const _convSIspeed = Dict("m/s" => 1.0, 
                          "kts" => kts_to_mps,
                          "km/hr"=>1000.0/3600.0,
                          "ft/s"=>ft_to_m/1.0)

function convertSpeed(value::Float64, units_in::AbstractString="m/s", units_out="m/s")
    dict = _convSIspeed
    return value * dict[units_in]/dict[units_out]
end

#Power
const _convSIpower = Dict("W"=>1.0,
                        "kW"=>1000.0,
                        "MW"=>1e6,
                        "hp"=>hp_to_W)

function convertPower(value::Float64, units_in::AbstractString="W", units_out="W")
    dict = _convSIpower
    return value * dict[units_in]/dict[units_out]
end

#Time
const _convSItime = Dict("s"=>1.0,
                        "h"=>3600.0,
                        "days"=>86400.0)

function convertTime(value::Float64, units_in::AbstractString="s", units_out="s")
    dict = _convSItime
    return value * dict[units_in]/dict[units_out]
end

#Angle
const _convSIangle = Dict("rad" => 1.0, "deg"=> deg_to_rad)

function convertAngle(value::Float64, units_in::AbstractString="rad", units_out="rad")
    dict = _convSIangle
    return value * dict[units_in]/dict[units_out]
end

function convertTemp(value::Float64, units_in::AbstractString="K", units_out = "K")
    units_in = uppercase(units_in)
    units_out = uppercase(units_out)
    if units_in == units_out
        return value
    else
        # First just convert to K
        if units_in == "K"
            temp = value
        elseif units_in == "C"
            temp =  value + 273.15
        elseif units_in == "F"
            temp = (value - 32.0) * 100.0/180.0 + 273.5
        elseif units_in == "R"
            temp = value /1.8
        end
        # Convert to desired output unit
        if units_out == "K"
            return temp
        elseif units_out == "C"
            return temp - 273.15
        elseif units_out == "F"
            return (temp - 273.15) * 180.0/100.0 + 32
        elseif units_out == "R"
            return temp*1.8
        end
    end

end

_allowed_units = vcat(([keys(x) for x in [_convSIdist,
                                          _convSImass,
                                          _convSIforce,
                                          _convSIspeed,
                                          _convSIpower,
                                          _convSItime,
                                          _convSIangle]]...)...)



