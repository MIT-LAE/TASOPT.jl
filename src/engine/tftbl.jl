"""
`tftbl.jl` loads in the E3 map data for bilinear interpolation in tfmap.jl
"""
struct compressorTbl
    # ---- Design values
    Nb_des::AbstractFloat
    Rline_des::AbstractFloat # Currently unused
    mb_des::AbstractFloat
    pr_des::AbstractFloat
    eff_is_des::AbstractFloat # Isentropic design efficiency

    # ---- Interpolators (Nominal and Gradient)
    Nc_itp::natNeigh.NaturalNeighboursInterpolant  
    Nc_grad::natNeigh.NaturalNeighboursDifferentiator    
    eff_itp::natNeigh.NaturalNeighboursInterpolant  
    eff_grad::natNeigh.NaturalNeighboursDifferentiator    
end 


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
function create_map_struct(map_name::String; method::natNeigh.AbstractInterpolator=natNeigh.Sibson(), Nres::Int=50, Rres::Int=50)
    # Load map data from TOML file
    println("Generatring Map: " * map_name)
    map = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/",map_name))

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

    # Create the gradient calculators for each output variable
    println("Creating final gradients")
    N_gradients::natNeigh.NaturalNeighboursDifferentiator = natNeigh.differentiate(N_interpolator, 1)
    E_gradients::natNeigh.NaturalNeighboursDifferentiator = natNeigh.differentiate(E_interpolator, 1)

    # Get the design quantities for each relevant input/output
    Nb_des::Float64 = map["design"]["des_Nc"]
    Rline_des::Float64 = map["design"]["des_Rline"]
    
    N_ind::Int = findall(x->isapprox(x,Nb_des), map["grid"]["grid_Nc"])[1]
    R_ind::Int = findall(x->isapprox(x,Rline_des), map["grid"]["grid_Rline"])[1]

    mb_des::Float64  = map["characteristic"]["massflow"][N_ind][R_ind]
    pr_des::Float64  = map["characteristic"]["pressure_ratio"][N_ind][R_ind]
    eff_des::Float64 = map["characteristic"]["isen_efficiency"][N_ind][R_ind]

    # Construct the compressorTbl struct for the given map
    mapStruct = compressorTbl(
        Nb_des,
        Rline_des,
        mb_des,
        pr_des,
        eff_des,
        N_interpolator,
        N_gradients,
        E_interpolator,
        E_gradients
    )

    println("Done Generating\n\n")
    return mapStruct
end # create_map_struct


if isfile(joinpath(__TASOPTroot__,"engine/data/E3fan_gridded.toml"))
    fanDict = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/E3fan_gridded.toml"))
else
    fanDict = create_map_struct("E3fan.toml"; method=natNeigh.Triangle(0; allow_cache=false), Nres=250, Rres=250)
end






const E3lpc = create_map_struct("E3lpc.toml"; method=natNeigh.Sibson(0), Nres=50, Rres=50)
const E3hpc = create_map_struct("E3hpc.toml"; method=natNeigh.Triangle(0; allow_cache=false), Nres=25, Rres=25)


