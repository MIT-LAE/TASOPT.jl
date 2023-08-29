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



