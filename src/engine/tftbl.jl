"""
`tftbl.jl` loads in the E3 map data for bilinear interpolation in tfmap.jl
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


function gen2DGrid(itp, xLow, xHigh, yLow, yHigh, Nx, Ny; method=natNeigh.Sibson(1))
    xg = range(xLow, xHigh, Nx)
    yg = range(yLow, yHigh, Ny)

    dx = (xHigh - xLow) / (Nx-1)
    dy = (yHigh - yLow) / (Ny-1)

    xVec = vec([i for i in xg, _ in yg])
    yVec = vec([j for _ in xg, j in yg])

    fVec = itp(xVec, yVec; method=method)

    fGridded = reshape(fVec, (length(xg), length(yg)))

    return fGridded, xg, yg, dx, dy
end # gen2DGrid


function gen2DDerivatives(grid, dx, dy, Nx, Ny)
    xDeriv = stack(fill(0., Nx, Ny))
    yDeriv = stack(fill(0., Nx, Ny))

    for i in 1:Nx
        for j in 1:Ny
            if i == 1
                xDeriv[i,j] = -3. * grid[i,j] + 4. * grid[i+1, j] - grid[i+2, j]
            elseif i == Nx
                xDeriv[i,j] = 3. * grid[i,j] - 4. * grid[i-1, j] + grid[i-2, j]
            else
                xDeriv[i,j] = grid[i+1,j] - grid[i-1,j]
            end

            if j == 1
                yDeriv[i,j] = -3. * grid[i,j] + 4. * grid[i, j+1] - grid[i, j+2]
            elseif j == Ny
                yDeriv[i,j] =  3. * grid[i,j] - 4. * grid[i, j-1] + grid[i, j-2]
            else
                yDeriv[i,j] = grid[i,j+1] - grid[i,j-1]
            end
        end
    end

    xDeriv = xDeriv ./ (2*dx)
    yDeriv = yDeriv ./ (2*dy)

    return xDeriv, yDeriv
end # gen2DDerivatives


"""
    create_map_struct(map_name; method=natNeigh.Sibson())

Creates a struct representative of a compressor map given toml inputs like those included
in the data folder (see E3fan.toml for an example). The struct includes design values of
corrected speed, Rline, mass flow, pressure ratio, and isentropic efficiency. It also 
includes four functions meant to act as interpolators for the nominal values and gradients
for correected speed and efficiency.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `map_name`:  compressor map data file name in format name.toml
    - `method`:    interpolant method; defaults to Sibson from the NaturalNeighbours package

    **Outputs:**
    - `mapStruct`: compressorTbl representation of the map

"""
function create_map_struct(map_name::String; method::natNeigh.AbstractInterpolator=natNeigh.Sibson(), Nres::Int=250, Rres::Int=250, Mres::Int=100, Pres::Int=100)
    # Load map data from TOML file
    println("Generatring Map: " * map_name)
    map = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/",map_name))
    hash_val = hash(map)

    # Need to specify each point in the N-R plane (i.e. all combinations of N_i and R_i)
    # To do so, create a vector with all x (N) components and a vector with all y (R) components
    Nvec::Vector{Float64} = vec([i for i in map["grid"]["grid_Nc"], _ in map["grid"]["grid_Rline"]])
    Rvec::Vector{Float64} = vec([j for _ in map["grid"]["grid_Nc"], j in map["grid"]["grid_Rline"]])

    # Load in the mass flow, pressure ratio, and isentropic efficiency values in the same order as the points are specified above
    Mvec::Vector{Float64} = vec([map["characteristic"]["massflow"][i][j]        for i in 1:length(map["grid"]["grid_Nc"]), j in 1:length(map["grid"]["grid_Rline"])])
    Pvec::Vector{Float64} = vec([map["characteristic"]["pressure_ratio"][i][j]  for i in 1:length(map["grid"]["grid_Nc"]), j in 1:length(map["grid"]["grid_Rline"])])
    Evec::Vector{Float64} = vec([map["characteristic"]["isen_efficiency"][i][j] for i in 1:length(map["grid"]["grid_Nc"]), j in 1:length(map["grid"]["grid_Rline"])])

    # Create Natural Neighbour interpolators for each parameter in terms of N and R
    println("Creating Grid Interpolators")
    Mitp::natNeigh.NaturalNeighboursInterpolant = natNeigh.interpolate(Nvec, Rvec, Mvec; derivatives=true)
    Pitp::natNeigh.NaturalNeighboursInterpolant = natNeigh.interpolate(Nvec, Rvec, Pvec; derivatives=true)
    Eitp::natNeigh.NaturalNeighboursInterpolant = natNeigh.interpolate(Nvec, Rvec, Evec; derivatives=true)
    
    # Create a finer grid of test points; convert them to interpolant input using the same method as Nvec and Rvec (above)
    N_testGrid::AbstractRange{Float64} = range(minimum(map["grid"]["grid_Nc"]),    maximum(map["grid"]["grid_Nc"]),    Nres)
    R_testGrid::AbstractRange{Float64} = range(minimum(map["grid"]["grid_Rline"]), maximum(map["grid"]["grid_Rline"]), Rres)

    N_testVec::Vector{Float64} = vec([i for i in N_testGrid, _ in R_testGrid])
    R_testVec::Vector{Float64} = vec([j for _ in N_testGrid, j in R_testGrid])

    # Get the mass flow, pressure ratio, and efficiency values at the finer grid using the interpolators
    println("Starting gridding interpolation")
    M_testVec::Vector{Float64} = Mitp(N_testVec, R_testVec; method=method, parallel=false)
    P_testVec::Vector{Float64} = Pitp(N_testVec, R_testVec; method=method, parallel=false)
    E_testVec::Vector{Float64} = Eitp(N_testVec, R_testVec; method=method, parallel=false)

    # Create direct interpolators from mb, PR -> Nb and efficiency
    println("Creating final interpolators")
    N_interpolator::natNeigh.NaturalNeighboursInterpolant = natNeigh.interpolate(M_testVec, P_testVec, N_testVec; derivatives=true)
    E_interpolator::natNeigh.NaturalNeighboursInterpolant = natNeigh.interpolate(M_testVec, P_testVec, E_testVec; derivatives=true)

    # Use the interpolators to create a regular grid of N and efficiency values
    println("Gridding data")
    N_gridded, Mg, Pg, dm, dp = gen2DGrid(N_interpolator, minimum(Mvec), maximum(Mvec), minimum(Pvec), maximum(Pvec), Mres, Pres; method=natNeigh.Sibson(1))
    E_gridded, _,  _,  _,  _  = gen2DGrid(E_interpolator, minimum(Mvec), maximum(Mvec), minimum(Pvec), maximum(Pvec), Mres, Pres; method=natNeigh.Sibson(1)) 

    # Calculate the partial derivatives of N and efficiency to mass flow and pressure ratio on the same grid
    println("Calculating partials")
    dNdm, dNdp = gen2DDerivatives(N_gridded, dm, dp, Mres, Pres)
    dEdm, dEdp = gen2DDerivatives(E_gridded, dm, dp, Mres, Pres)

    # Get the design quantities for each relevant input/output
    println("Getting design values")
    Nb_des::Float64 = map["design"]["des_Nc"]
    Rline_des::Float64 = map["design"]["des_Rline"]
    
    N_ind::Int = findall(x->isapprox(x,Nb_des), map["grid"]["grid_Nc"])[1]
    R_ind::Int = findall(x->isapprox(x,Rline_des), map["grid"]["grid_Rline"])[1]

    mb_des::Float64  = map["characteristic"]["massflow"][N_ind][R_ind]
    pr_des::Float64  = map["characteristic"]["pressure_ratio"][N_ind][R_ind]
    eff_des::Float64 = map["characteristic"]["isen_efficiency"][N_ind][R_ind]

    # Construct the compressorTbl struct for the given map
    println("Constructing compressorTbl")
    mapStruct = compressorTbl(
        Nb_des,
        Rline_des,
        mb_des,
        pr_des,
        eff_des,
        Mg,
        Pg,
        dm,
        dp,
        Mres,
        Pres,
        N_gridded,
        dNdm,
        dNdp,
        E_gridded,
        dEdm,
        dEdp,
        hash_val,
    )

    println("Done Generating\n\n")
    return mapStruct
end # create_map_struct


function checkMapHash(gridded_map_name, data_map_name)
    println("Checking hash for"*data_map_name)
    genComp = true
    compDict = Dict()
    cacheDict = Dict()
    if isfile(joinpath(__TASOPTroot__,"engine/data/",gridded_map_name))
        compDict = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/",data_map_name))
        cacheDict = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/",gridded_map_name))
        if hash(compDict) == cacheDict["hash_val"]
            genComp = false
        end
    end

    return genComp, cacheDict
end # checkMapHash


function map_startup(gridded_map_name, data_map_name, map_gridding_method; Nres=250, Rres=250, Mres=100, Pres=100)
    genComp, compDict = checkMapHash(gridded_map_name, data_map_name)

    if genComp
        # Generate new fan data and save to a file
        compStruct = create_map_struct(data_map_name; method=map_gridding_method, Nres=Nres, Rres=Rres, Mres=Mres, Pres=Pres)
        compTbl_to_TOML(compStruct, gridded_map_name)    
    else
        compStruct = dict_to_compTbl(compDict)
    end

    return compStruct
end # map_startup


# Create compressorTbl objects as constants to prevent re-generating all the maps every run
const E3fan = map_startup("E3fan_gridded.toml", "E3fan.toml", natNeigh.Triangle(0; allow_cache=false); Nres=300, Rres=300, Mres=100, Pres=100)
const E3lpc = map_startup("E3lpc_gridded.toml", "E3lpc.toml", natNeigh.Sibson(0);                      Nres=300, Rres=300, Mres=100, Pres=100)
const E3hpc = map_startup("E3hpc_gridded.toml", "E3hpc.toml", natNeigh.Triangle(0; allow_cache=false); Nres=300, Rres=300, Mres=100, Pres=100)
