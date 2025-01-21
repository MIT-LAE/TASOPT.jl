"""
    pointInMap(m, p, surge, windmill, bottom, top)

Checks whether a mass flow / pressure ratio map is in the domain of the map
defined by N and Rline. Assumes the boundary lists are pre-sorted

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `m`:        mass flow test value
    - `p`:        pressure ratio test value
    - `surge`:    matrix of m/p coordinates defining surge line
    - `windmill`: matrix of m/p coordinates defining windmill line
    - `bottom`:   matrix of m/p coordinates defining topmost N line
    - `top`:      matrix of m/p coordinates defining bottommost N line

    **Outputs:**
    - `inMap`:    boolean; true if point is in map; false otherwise

"""
function pointInMap(m, p, surge, windmill, bottom, top)
    # Find the transition point on the x axis between the upper and lower bounding lines of the map
    mMid_top = maximum(surge[:,1])
    mMid_bot = maximum(bottom[:,1])

    top_bound = surge
    bot_bound = bottom
    # println(m, " ", p)

    # Set the upper and lower bound appropriately for the given m value
    if m > mMid_top
        top_bound = top
        # println("top")
    else
        # println("surge")
    end

    if m > mMid_bot
        bot_bound = windmill
        # println("windmill")
    else
        # println("bottom")
    end

    # Get the indices of points the m value lies between for interpolation
    indUpp_top = searchsortedfirst(top_bound[:,1], m)
    indLow_top = max(indUpp_top - 1, 1)
    indUpp_top = min(indUpp_top, length(top_bound[:,1]))
    low_top = top_bound[indLow_top,:]
    upp_top = top_bound[indUpp_top,:]

    indUpp_bot = searchsortedfirst(bot_bound[:,1], m)
    indLow_bot = max(indUpp_bot - 1, 1)
    indUpp_bot = min(indUpp_bot, length(bot_bound[:,1]))
    # println(indLow_bot + 1, " ", length(bot_bound[:,1]), " ", indUpp_bot)
    low_bot = bot_bound[indLow_bot,:]
    upp_bot = bot_bound[indUpp_bot,:]

    # Calculate the minimum/maximum pr value required to be in map using linear interpolation
    if indLow_top != indUpp_top
        dm_top = upp_top[1] - low_top[1]
        wLow_top = (upp_top[1] - m) / dm_top
        wUpp_top = (m - low_top[1]) / dm_top
        pmax = low_top[2] * wLow_top + upp_top[2] * wUpp_top
    else
        pmax = upp_top[2]
    end

    if indLow_bot != indUpp_bot
        dm_bot = upp_bot[1] - low_bot[1]
        wLow_bot = (upp_bot[1] - m) / dm_bot
        wUpp_bot = (m - low_bot[1]) / dm_bot
        pmin = low_bot[2] * wLow_bot + upp_bot[2] * wUpp_bot
    else
        pmin = upp_bot[2]
    end

    # If greater than pmax or lower than pmin, the point isn't in the map
    if p < pmin || p > pmax
        return false
    else
        return true
    end
    
end # pointInMap


"""
    gen2DGrid(itp, xLow, xHigh, yLow, yHigh, Nx, Ny; method=natNeigh.Sibson(1))

Creates a struct representative of a compressor map given toml inputs like those included
in the data folder (see E3fan.toml for an example). The struct includes design values of
corrected speed, Rline, mass flow, pressure ratio, and isentropic efficiency. It also 
includes four functions meant to act as interpolators for the nominal values and gradients
for correected speed and efficiency.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `itp`: NaturalNeighbours.NaturalNeighboursInterpolant object such that (x, y) -> f, where f is the interpolated data output
    - `xLow`: minimum value of the x grid
    - `xHigh`: maximum value of the x grid
    - `yLow`: minimum value of the y grid
    - `yHigh`: maximum value of the y grid
    - `Nx`: number of x grid points
    - `Ny`: number of y grid points
    - `method`: interpolant method; defaults to Sibson-1 from the NaturalNeighbours package

    **Outputs:**
    - `fGridded`: Nx by Ny matrix of regularly gridded function values
    - `xg`: vector of x grid values
    - `yg`: vector of y grid values
    - `dx`: x grid spacing
    - `dy`: y grid spacing

"""
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


"""
    gen2DDerivatives(grid, dx, dy, Nx, Ny)

Uses central finite difference method to estimate x and y derivatives of grid values. At boundaries, uses second-order accurate forward/backwards derivatives.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `itp`: NaturalNeighbours.NaturalNeighboursInterpolant object such that (x, y) -> f, where f is the interpolated data output
    - `dx`: x grid spacing
    - `dy`: y grid spacing
    - `Nx`: number of x grid points
    - `Ny`: number of y grid points

    **Outputs:**
    - `xDeriv`: Nx by Ny matrix of x derivatives
    - `yDeriv`: Nx by Ny matrix of y derivatives

"""
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
function create_map_struct(map_name::String, oob_consts::Vector{Float64}; method::natNeigh.AbstractInterpolator=natNeigh.Sibson(), Nres::Int=250, Rres::Int=250, Mres::Int=100, Pres::Int=100)
    # Load map data from TOML file
    println("Generatring Map: " * map_name)
    map = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/",map_name))
    hash_val = hash(map)

    # Get the design quantities for each relevant input/output
    println("Getting design values")
    Nb_des::Float64 = map["design"]["des_Nc"]
    Rline_des::Float64 = map["design"]["des_Rline"]
    
    N_ind::Int = findall(x->isapprox(x,Nb_des), map["grid"]["grid_Nc"])[1]
    R_ind::Int = findall(x->isapprox(x,Rline_des), map["grid"]["grid_Rline"])[1]

    mb_des::Float64  = map["characteristic"]["massflow"][N_ind][R_ind]
    pr_des::Float64  = map["characteristic"]["pressure_ratio"][N_ind][R_ind]
    eff_des::Float64 = map["characteristic"]["isen_efficiency"][N_ind][R_ind]

    # Get map boundaries in the mb - pr space
    surgex = stack(map["characteristic"]["massflow"])[1,:]
    surgey = stack(map["characteristic"]["pressure_ratio"])[1,:]
    surge = hcat(surgex, surgey)

    windmillx = stack(map["characteristic"]["massflow"])[end,:]
    windmilly = stack(map["characteristic"]["pressure_ratio"])[end,:]
    windmill = hcat(windmillx, windmilly)

    bottomx = stack(map["characteristic"]["massflow"])[:,1]
    bottomy = stack(map["characteristic"]["pressure_ratio"])[:,1]
    bottom = hcat(bottomx, bottomy)

    topx = stack(map["characteristic"]["massflow"])[:,end]
    topy = stack(map["characteristic"]["pressure_ratio"])[:,end]
    top = hcat(topx, topy)

    # Need to specify each point in the N-R plane (i.e. all combinations of N_i and R_i)
    # To do so, create a vector with all x (N) components and a vector with all y (R) components
    Nvec::Vector{Float64} = vec([i for i in map["grid"]["grid_Nc"], _ in map["grid"]["grid_Rline"]])
    Rvec::Vector{Float64} = vec([j for _ in map["grid"]["grid_Nc"], j in map["grid"]["grid_Rline"]])

    # Load in the mass flow, pressure ratio, and isentropic efficiency values in the same order as the points are specified above
    # Make sure to normalize the mass flow and pressure ratio by the design values to get them on a similar scale
    isen_to_poly = (pr, isen, g) -> (g-1)/g * log(pr) / log( (pr^((g-1)/g) - 1) / isen + 1 )
    Mvec::Vector{Float64} = vec([map["characteristic"]["massflow"][i][j]/mb_des        for i in 1:length(map["grid"]["grid_Nc"]), j in 1:length(map["grid"]["grid_Rline"])])
    Pvec::Vector{Float64} = vec([map["characteristic"]["pressure_ratio"][i][j]/pr_des  for i in 1:length(map["grid"]["grid_Nc"]), j in 1:length(map["grid"]["grid_Rline"])])
    
    # Convert from isentropic to polytropic efficiency here; gamma is considered constant at 1.4
    Evec::Vector{Float64} = vec([isen_to_poly(map["characteristic"]["pressure_ratio"][i][j], map["characteristic"]["isen_efficiency"][i][j], 1.4) for i in 1:length(map["grid"]["grid_Nc"]), j in 1:length(map["grid"]["grid_Rline"])])

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

    # Scale the mass flow and pressure back up post-interpolation
    Mg = Mg .* mb_des
    dm = dm  * mb_des

    Pg = Pg .* pr_des
    dp = dp  * pr_des

    # Replace out-of-map points with the Drela approximation
    for i in eachindex(Mg)
        for j in eachindex(Pg)
            if !pointInMap(Mg[i], Pg[j], surge, windmill, bottom, top)
                N_gridded[i,j], _, _ = Ncmap(Pg[j], Mg[i], pr_des, mb_des, Nb_des, oob_consts)
                E_gridded[i,j], _, _ = ecmap(Pg[j], Mg[i], pr_des, mb_des, oob_consts, isen_to_poly(pr_des, eff_des, 1.4), 0.0, 0.0)
            end
        end
    end

    # Calculate the partial derivatives of N and efficiency to mass flow and pressure ratio on the same grid
    println("Calculating partials")
    dNdm, dNdp = gen2DDerivatives(N_gridded, dm, dp, Mres, Pres)
    dEdm, dEdp = gen2DDerivatives(E_gridded, dm, dp, Mres, Pres)

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
        surge,
        windmill,
        bottom,
        top,
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


"""
    checkMapHash(gridded_map_name, data_map_name)

Loads gridded and NPSS map data from the `src/engine/data/` folder and checks whether the hash value stored in the gridded file is equal to the hash value of the
NPSS map data.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gridded_map_name`: name of the gridded data file
    - `data_map_name`: name of the NPSS map data file

    **Outputs:**
    - `genComp`: boolean; false if the hash matches, true if the hash does not
    - `cacheDict`: dictionary representation of the gridded data

"""
function checkMapHash(gridded_map_name, data_map_name)
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


"""
    map_startup(gridded_map_name, data_map_name, map_gridding_method; Nres=250, Rres=250, Mres=100, Pres=100)

Map loading and generating function, runs on the time of precompilation. Determines whether a gridded version of the map data already exists and loads it if so. Otherwise
the map is generated from the NPSS data using `create_map_struct`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gridded_map_name`: `.toml` file name of the gridded data
    - `data_map_name`: `.toml` file name of the NPSS map data
    - `map_gridding_method`: NaturalNeighbours.NaturalNeighboursInterpolant to be used for map gridding (if needed)
    - `Nres`: resolution of corrected speed sampling; default 200
    - `Rres`: resolution of R-line sampling; default 200
    - `Mres`: resolution of corrected mass flow sampling; default 200
    - `Pres`: resolution of pressure ratio sampling; default 200

    **Outputs:**
    - `compStruct`: compressorTbl structure representing the loaded map

"""
function map_startup(gridded_map_name, data_map_name, map_gridding_method, oob_consts; Nres=200, Rres=200, Mres=200, Pres=200)
    genComp, compDict = checkMapHash(gridded_map_name, data_map_name)

    if genComp
        # Generate new fan data and save to a file
        compStruct = create_map_struct(data_map_name, oob_consts; method=map_gridding_method, Nres=Nres, Rres=Rres, Mres=Mres, Pres=Pres)
        compTbl_to_TOML(compStruct, gridded_map_name)    
    else
        compStruct = dict_to_compTbl(compDict)
    end

    return compStruct
end # map_startup


# Specify desired griding fineness; possibly make this a startup variable?
N, R, M, P = 200, 200, 400, 400


# Load Drela constants for out-of-bounds calculations
#                  a     b     k     mo     da    c    d    C     D
const Cmapf_oob = [3.50, 0.80, 0.03, 0.95, -0.50, 3.0, 6.0,  2.5, 15.0]
const Cmapl_oob = [1.90, 1.00, 0.03, 0.95, -0.20, 3.0, 5.5, 15.0,  1.0]
const Cmaph_oob = [1.75, 2.00, 0.03, 0.95, -0.35, 3.0, 5.0, 15.0,  1.0]

# Create compressorTbl objects as constants to prevent re-generating all the maps every run
const E3fan = map_startup("E3fan_gridded.toml", "E3fan.toml", natNeigh.Sibson(1), Cmapf_oob; Nres=N, Rres=R, Mres=M, Pres=P)
const E3lpc = map_startup("E3lpc_gridded.toml", "E3lpc.toml", natNeigh.Sibson(1), Cmapl_oob; Nres=N, Rres=R, Mres=M, Pres=P)
const E3hpc = map_startup("E3hpc_gridded.toml", "E3hpc.toml", natNeigh.Sibson(1), Cmaph_oob; Nres=N, Rres=R, Mres=M, Pres=P)

# Turbine map constants
const Tmapl = [0.15, 0.15]
const Tmaph = [0.15, 0.15]