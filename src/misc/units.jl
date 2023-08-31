#Dist
const _convSIdist = Dict("m" => 1.0, 
                        "km"=>1000.0, 
                        "nmi" => nmi_to_m,
                        "ft"=> ft_to_m, 
                        "in"=> in_to_m)

function convertDist(value::Float64, units_in::String, units_out="m")
    dict = _convSIdist
    return value * dict[units_in]/dict[units_out]
end

"""
    parse_unit_string(str::String)

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
function parse_unit_string(str::String)
    m = match(r"(\+?\-?\d+\.?e?\d*)\s*(\w*\/?\w*)", str)
    mag = parse(Float64, m[1])
    units = m[2]
    return mag, units
end

function get_unit_dim(str::String)
    m = match(r"([^\d]+)(\d*)", str)
    unit = m[1]
    if m[2] == ""
        dim = 1.0
    else
        dim = parse(Float64, m[2])
    end
    return unit, dim
end

function convertArea(value::Float64, units_in::String, units_out="m2")
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

function convertVolume(value::Float64, units_in::String, units_out="m3")
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
                         "lbm"=>1.0/2.205)

function convertMass(value::Float64, units_in::String, units_out="kg")
    dict = _convSImass
    return value * dict[units_in]/dict[units_out]
end

#Force
const _convSIforce = Dict("N" => 1.0,
                          "kN" => 1000.0,
                          "lbf" => lbf_to_N)       

function convertForce(value::Float64, units_in::String, units_out="N")
    dict = _convSIforce
    return value * dict[units_in]/dict[units_out]
end

#Speed
const _convSIspeed = Dict("m/s" => 1.0, 
                          "kts" => kts_to_mps,
                          "km/hr"=>1000.0/3600.0,
                          "ft/s"=>ft_to_m/1.0)

function convertSpeed(value::Float64, units_in::String, units_out="m/s")
    dict = _convSIspeed
    return value * dict[units_in]/dict[units_out]
end

#Power
const _convSIpower = Dict("W"=>1.0,
                        "kW"=>1000.0,
                        "MW"=>1e6,
                        "hp"=>hp_to_W)

function convertPower(value::Float64, units_in::String, units_out="W")
    dict = _convSIpower
    return value * dict[units_in]/dict[units_out]
end

#Angle
const _convSIangle = Dict("rad" => 1.0, "deg"=> deg_to_rad)

function convertAngle(value::Float64, units_in::String, units_out="deg")
    dict = _convSIangle
    return value * dict[units_in]/dict[units_out]
end

_allowed_units = vcat(([keys(x) for x in [_convSIdist,
                                          _convSImass,
                                          _convSIforce,
                                          _convSIspeed,
                                          _convSIpower,
                                          _convSIangle]]...)...)



