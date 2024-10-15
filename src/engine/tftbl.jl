"""
`tftbl.jl` loads in the E3 map data for bilinear interpolation in tfmap.jl
"""
# using TOML
# import ..TASOPT: __TASOPTroot__

struct compressorTbl
    # ---- Independent variables
    mb_map::AbstractArray{AbstractFloat}
    pr_map::AbstractArray{AbstractFloat}

    # ---- Design values
    Nb_des::AbstractFloat
    Rline_des::AbstractFloat # Currently unused
    mb_des::AbstractFloat
    pr_des::AbstractFloat
    eff_is_des::AbstractFloat # Isentropic design efficiency

    # ---- Lookup tables
    Nb_tbl::AbstractArray{AbstractFloat}
    Rline_tbl::AbstractArray{AbstractFloat}
    eff_is_tbl::AbstractArray{AbstractFloat}
    eff_poly_tbl::AbstractArray{AbstractFloat}
end 

E3fan_toml = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/E3fan_tbl.toml"))
E3lpc_toml = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/E3lpc_tbl.toml"))
E3hpc_toml = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/E3hpc_tbl.toml"))

# E3fan_toml = TOML.parsefile("$(@__DIR__)/data/E3fan_tbl.toml")
# E3lpc_toml = TOML.parsefile("$(@__DIR__)/data/E3lpc_tbl.toml")
# E3hpc_toml = TOML.parsefile("$(@__DIR__)/data/E3hpc_tbl.toml")

const E3fan = compressorTbl(
    E3fan_toml["input"]["mb"],
    E3fan_toml["input"]["pr"],
    E3fan_toml["design"]["Nb_des"],
    E3fan_toml["design"]["Rline_des"],
    E3fan_toml["design"]["mb_des"],
    E3fan_toml["design"]["pr_des"],
    E3fan_toml["design"]["eff_des"],
    stack(E3fan_toml["output"]["Nb"]),
    stack(E3fan_toml["output"]["R"]),
    stack(E3fan_toml["output"]["Ei"]),
    stack(E3fan_toml["output"]["Ep"]),
)

const E3lpc = compressorTbl(
    E3lpc_toml["input"]["mb"],
    E3lpc_toml["input"]["pr"],
    E3lpc_toml["design"]["Nb_des"],
    E3lpc_toml["design"]["Rline_des"],
    E3lpc_toml["design"]["mb_des"],
    E3lpc_toml["design"]["pr_des"],
    E3lpc_toml["design"]["eff_des"],
    stack(E3lpc_toml["output"]["Nb"]),
    stack(E3lpc_toml["output"]["R"]),
    stack(E3lpc_toml["output"]["Ei"]),
    stack(E3lpc_toml["output"]["Ep"]),
)

const E3hpc = compressorTbl(
    E3hpc_toml["input"]["mb"],
    E3hpc_toml["input"]["pr"],
    E3hpc_toml["design"]["Nb_des"],
    E3hpc_toml["design"]["Rline_des"],
    E3hpc_toml["design"]["mb_des"],
    E3hpc_toml["design"]["pr_des"],
    E3hpc_toml["design"]["eff_des"],
    stack(E3hpc_toml["output"]["Nb"]),
    stack(E3hpc_toml["output"]["R"]),
    stack(E3hpc_toml["output"]["Ei"]),
    stack(E3hpc_toml["output"]["Ep"]),
)