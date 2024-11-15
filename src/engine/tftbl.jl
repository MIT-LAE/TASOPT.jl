"""
`tftbl.jl` loads in the E3 map data for bilinear interpolation in tfmap.jl
"""
# using TOML
# import ..TASOPT: __TASOPTroot__

# struct compressorTbl
#     # ---- Independent variables
#     mb_map::AbstractArray{AbstractFloat}
#     pr_map::AbstractArray{AbstractFloat}

#     # ---- Design values
#     Nb_des::AbstractFloat
#     Rline_des::AbstractFloat # Currently unused
#     mb_des::AbstractFloat
#     pr_des::AbstractFloat
#     eff_is_des::AbstractFloat # Isentropic design efficiency

#     # ---- Lookup tables
#     Nb_tbl::AbstractArray{AbstractFloat}
#     Rline_tbl::AbstractArray{AbstractFloat}
#     eff_is_tbl::AbstractArray{AbstractFloat}
#     eff_poly_tbl::AbstractArray{AbstractFloat}
# end 

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
function create_map_struct(map_name; method=natNeigh.Sibson(), Nres=50, Rres=50)
    # Load map data from TOML file
    println("Generatring Map: " * map_name)
    map = TOML.parsefile(joinpath(__TASOPTroot__,"engine/data/",map_name))

    # Need to specify each point in the N-R plane (i.e. all combinations of N_i and R_i)
    # To do so, create a vector with all x (N) components and a vector with all y (R) components
    Nvec = vec([i for i in map["grid"]["grid_Nc"], _ in map["grid"]["grid_Rline"]])
    Rvec = vec([j for _ in map["grid"]["grid_Nc"], j in map["grid"]["grid_Rline"]])

    # Load in the mass flow, pressure ratio, and isentropic efficiency values in the same order as the points are specified above
    Mvec = vec([map["characteristic"]["massflow"][i][j]        for i in 1:length(map["grid"]["grid_Nc"]), j in 1:length(map["grid"]["grid_Rline"])])
    Pvec = vec([map["characteristic"]["pressure_ratio"][i][j]  for i in 1:length(map["grid"]["grid_Nc"]), j in 1:length(map["grid"]["grid_Rline"])])
    Evec = vec([map["characteristic"]["isen_efficiency"][i][j] for i in 1:length(map["grid"]["grid_Nc"]), j in 1:length(map["grid"]["grid_Rline"])])

    # Create Natural Neighbour interpolators for each parameter in terms of N and R
    println("Creating Grid Interpolators")
    Mitp = natNeigh.interpolate(Nvec, Rvec, Mvec; derivatives=true)
    Pitp = natNeigh.interpolate(Nvec, Rvec, Pvec; derivatives=true)
    Eitp = natNeigh.interpolate(Nvec, Rvec, Evec; derivatives=true)
    
    # Create a finer grid of test points; convert them to interpolant input using the same method as Nvec and Rvec (above)
    N_testGrid = range(minimum(map["grid"]["grid_Nc"]),    maximum(map["grid"]["grid_Nc"]),    Nres)
    R_testGrid = range(minimum(map["grid"]["grid_Rline"]), maximum(map["grid"]["grid_Rline"]), Rres)

    N_testVec = vec([i for i in N_testGrid, _ in R_testGrid])
    R_testVec = vec([j for _ in N_testGrid, j in R_testGrid])

    # Get the mass flow, pressure ratio, and efficiency values at the finer grid using the interpolators
    println("Starting gridding interpolation")
    M_testVec = Mitp(N_testVec, R_testVec; method=method, parallel=false)
    P_testVec = Pitp(N_testVec, R_testVec; method=method, parallel=false)
    E_testVec = Eitp(N_testVec, R_testVec; method=method, parallel=false)

    # Create direct interpolators from mb, PR -> Nb and efficiency
    println("Creating final interpolators")
    N_interpolator = natNeigh.interpolate(M_testVec, P_testVec, N_testVec; derivatives=true)
    E_interpolator = natNeigh.interpolate(M_testVec, P_testVec, E_testVec; derivatives=true)

    # Create the gradient calculators for each output variable
    println("Creating final gradients")
    N_gradients = natNeigh.differentiate(N_interpolator, 1)
    E_gradients = natNeigh.differentiate(E_interpolator, 1)

    # Get the design quantities for each relevant input/output
    Nb_des = map["design"]["des_Nc"]
    Rline_des = map["design"]["des_Rline"]
    
    N_ind = findall(x->isapprox(x,Nb_des), map["grid"]["grid_Nc"])[1]
    R_ind = findall(x->isapprox(x,Rline_des), map["grid"]["grid_Rline"])[1]

    mb_des  = map["characteristic"]["massflow"][N_ind][R_ind]
    pr_des  = map["characteristic"]["pressure_ratio"][N_ind][R_ind]
    eff_des = map["characteristic"]["isen_efficiency"][N_ind][R_ind]

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


const E3fan = create_map_struct("E3fan.toml"; method=natNeigh.Triangle(; allow_cache=false), Nres=50, Rres=50)
const E3lpc = create_map_struct("E3lpc.toml"; method=natNeigh.Sibson(), Nres=50, Rres=50)
const E3hpc = create_map_struct("E3hpc.toml"; method=natNeigh.Triangle(; allow_cache=false), Nres=25, Rres=25)


function evalNc(map, m, p; method=natNeigh.Sibson(1))
    return map.Nc_itp(m, p; method=method)
end

function evalEff(map, m, p; method=natNeigh.Sibson(1))
    return map.eff_itp(m, p; method=method)
end

function evalNcGrad(map, m, p)
    return map.Nc_grad(m, p)
end

function evalEffGrad(map, m, p)
    return map.eff_grad(m, p)
end