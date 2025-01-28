"""
`tftbl.jl` loads in the E3 map data for bilinear interpolation in tfmap.jl
"""

"""
    compressorTbl

A type representing a compressor map, parameterized in terms of mass flow and pressure ratio. Data is either generated from NPSS data using `create_map_struct` or is loaded
directly. In either case, the data files are stored in `src/engine/data`.

# Fields:
- `Nb_des::Float64`: Design corrected speed
- `Rline_des::Float64`: Design R-line value; currently unused
- `mb_des::Float64`: Design corrected mass flow
- `pr_des::Float64`: Design pressure ratio
- `eff_is_des::Float64`: Design isentropic efficiency

- `mbGrid::Vector{Float64}`: Grid of mass flow points at which output parameters are defined
- `prGrid::Vector{Float64}`: Grid of pressure ratio points at which output parameters are defined
- `dm::Float64`: Mass flow grid spacing
- `dp::Float64`: Pressure ratio grid spacing
- `Nm::Int`: Number of mass flow grid points
- `Np::Int`: Number of pressure ratio grid points

- `surge::Matrix{Float64}`:    matrix of m/p coordinates defining surge line
- `windmill::Matrix{Float64}`: matrix of m/p coordinates defining windmill line
- `bottom::Matrix{Float64}`:   matrix of m/p coordinates defining topmost N line
- `top::Matrix{Float64}`:      matrix of m/p coordinates defining bottommost N line

- `mask::Matrix{Bool}`

- `Nb_nom::Matrix{Float64}`: Nm by Np matrix of corrected speed values
- `Nb_mb::Matrix{Float64}`: Nm by Np matrix of corrected speed derivatives with respect to mass flow
- `Nb_pr::Matrix{Float64}`: Nm by Np matrix of corrected speed derivatives with respect to pressure ratio

- `Eff_nom::Matrix{Float64}`: Nm by Np matrix of polytropic efficiencies
- `Eff_mb::Matrix{Float64}`: Nm by Np matrix of polytropic efficiency with respect to mass flow
- `Eff_pr::Matrix{Float64}`: Nm by Np matrix of polytropic efficiency with respect to pressure ratio

- `hash_val::Integer`: Hash value corresponding to NPSS map data; tracks whether map data has been updated
"""
struct compressorTbl
    # ---- Design values
    Nb_des::Float64
    Rline_des::Float64 # Currently unused
    mb_des::Float64
    pr_des::Float64
    eff_is_des::Float64 # Isentropic design efficiency

    # ---- mb and PR grids and fineness variables
    mbGrid::Vector{Float64}
    prGrid::Vector{Float64}
    dm::Float64
    dp::Float64
    Nm::Int
    Np::Int

    # ---- Map boundaries
    surge::Matrix{Float64}
    windmill::Matrix{Float64}
    bottom::Matrix{Float64}
    top::Matrix{Float64}

    # ---- Mask of in-map points (1=in; 0=out)
    mask::Matrix{Bool}

    # ---- Nc grids
    Nb_nom::Matrix{Float64}
    Nb_mb::Matrix{Float64}
    Nb_pr::Matrix{Float64}

    # ---- Eff grids
    Eff_nom::Matrix{Float64}
    Eff_mb::Matrix{Float64}
    Eff_pr::Matrix{Float64}

    # ---- Hash
    hash_val::Integer
end 


"""
    compTbl_to_TOML(compressor::compressorTbl, saveName)

Converts a `compressorTbl` object to a TOML-saveable format and saves it as a TOML file to the `src/engine/data/` folder as `saveName`

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `compressor::compressorTbl`:  compressorTbl structure
    - `saveName`:    string file name for the saved data; should probably end in `.toml`
"""
function compTbl_to_TOML(compressor::compressorTbl, saveName)
    # Convert the struct fields to a dictionary
    compressor_dict = Dict(
        :Nb_des => compressor.Nb_des,
        :Rline_des => compressor.Rline_des,
        :mb_des => compressor.mb_des,
        :pr_des => compressor.pr_des,
        :eff_is_des => compressor.eff_is_des,
        :mbGrid => compressor.mbGrid,
        :prGrid => compressor.prGrid,
        :dm => compressor.dm,
        :dp => compressor.dp,
        :Nm => compressor.Nm,
        :Np => compressor.Np,
        :surge => [compressor.surge[i, :] for i in 1:size(compressor.surge, 1)],
        :windmill => [compressor.windmill[i, :] for i in 1:size(compressor.windmill, 1)],
        :bottom => [compressor.bottom[i, :] for i in 1:size(compressor.bottom, 1)],
        :top => [compressor.top[i, :] for i in 1:size(compressor.top, 1)],
        :mask => [compressor.mask[i, :] for i in 1:size(compressor.mask, 1)],
        :Nb_nom => [compressor.Nb_nom[i, :] for i in 1:size(compressor.Nb_nom, 1)],
        :Nb_mb => [compressor.Nb_mb[i, :] for i in 1:size(compressor.Nb_mb, 1)],
        :Nb_pr => [compressor.Nb_pr[i, :] for i in 1:size(compressor.Nb_pr, 1)],
        :Eff_nom => [compressor.Eff_nom[i, :] for i in 1:size(compressor.Eff_nom, 1)],
        :Eff_mb => [compressor.Eff_mb[i, :] for i in 1:size(compressor.Eff_mb, 1)],
        :Eff_pr => [compressor.Eff_pr[i, :] for i in 1:size(compressor.Eff_pr, 1)],
        :hash_val => compressor.hash_val
    )

    open(joinpath(__TASOPTroot__,"engine/data/",saveName), "w") do io
        TOML.print(io, compressor_dict)
    end
end # compTbl_to_TOML


"""
    dict_to_compTbl(compressor_dict)

Converts a compressorTbl dictionary (created when reading in a XXXX_gridded.toml file) to a compressorTbl object

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `compressor_dict`:  dictionary containing the compressor_dict data

    **Outputs:**
    - `compressor_struct`: compressorTbl representation of the data
"""
function dict_to_compTbl(compressor_dict)
    compressor_struct = compressorTbl(
        compressor_dict["Nb_des"],
        compressor_dict["Rline_des"],
        compressor_dict["mb_des"],
        compressor_dict["pr_des"],
        compressor_dict["eff_is_des"],
        compressor_dict["mbGrid"],
        compressor_dict["prGrid"],
        compressor_dict["dm"],
        compressor_dict["dp"],
        compressor_dict["Nm"],
        compressor_dict["Np"],
        transpose(hcat(compressor_dict["surge"]...)),
        transpose(hcat(compressor_dict["windmill"]...)),
        transpose(hcat(compressor_dict["bottom"]...)),
        transpose(hcat(compressor_dict["top"]...)),
        transpose(hcat(compressor_dict["mask"]...)),
        transpose(hcat(compressor_dict["Nb_nom"]...)),
        transpose(hcat(compressor_dict["Nb_mb"]...)),
        transpose(hcat(compressor_dict["Nb_pr"]...)),
        transpose(hcat(compressor_dict["Eff_nom"]...)),
        transpose(hcat(compressor_dict["Eff_mb"]...)),
        transpose(hcat(compressor_dict["Eff_pr"]...)),
        compressor_dict["hash_val"]
    )
    
    return compressor_struct
end # dict_to_compTbl