"""
    airfoil

A type representing a database of pre-computed airfoil data for a single Reyolds number and a range of Mach numbers, sectional lift coefficients, and thickness-to-chord ratios.
By default, this is the original TASOPT transonic airfoil, as modelled by M. Drela in MSES and stored in `src/air/`.

Overloads Base.show to print a summary of the `airfoil` model.

# Fields:
- `Ma::AbstractVector{Float64}` :  Mach nos.
- `cl::AbstractVector{Float64}` :  Sectional lift coefficients.
- `τ::AbstractVector{Float64}` :  Thickness-to-chord ratios.
- `Re::Float64` :  Reynolds number.

- `A::AbstractArray{Float64}`: Multi-dimensional array of aero performance data.

Various views of the data:
- `A_M::AbstractArray{Float64}`:
- `A_τ::AbstractArray{Float64}`:
- `A_cl::AbstractArray{Float64}`:
- `A_M_τ::AbstractArray{Float64}`:
- `A_M_cl::AbstractArray{Float64}`:
- `A_cl_τ::AbstractArray{Float64}`:
- `A_M_cl_τ::AbstractArray{Float64}`:

See also [`airfun`](@ref) and [`airtable`](@ref).
"""
struct airfoil
    Ma::AbstractVector{Float64}
    cl::AbstractVector{Float64}
    τ::AbstractVector{Float64}
    Re::Float64 # Data assumed for a single Re
  
    A::AbstractArray{Float64} # Airfoil aero data 
    
    A_M::AbstractArray{Float64}
    A_τ::AbstractArray{Float64}
    A_cl::AbstractArray{Float64}
    A_M_τ::AbstractArray{Float64}
    A_M_cl::AbstractArray{Float64}
    A_cl_τ::AbstractArray{Float64}
    A_M_cl_τ::AbstractArray{Float64}
end 

function Base.show(io::IO, airf::airfoil)
    printstyled(io, "Airfoil section database with:\n"; color=:bold)
    println(io, 
    "Re = ",airf.Re,"\n",
    "Ma = (", airf.Ma[1],", ", airf.Ma[end],")\n",
    "cl = (", airf.cl[1], ", ",airf.cl[end],")\n",
    "τ  = (", airf.τ[1], ", ",airf.τ[end],")\n"
    )
end

# using PythonPlot
# function plot(airf::airfoil)
#     fig, ax = plt.subplots(2,1, sharex = true)
#     ax[1].plot(airf.cl, airf.A[end, :, :, 1] + airf.A[end, :, :, 2], label = airf.τ)
#     ax[1].legend(title="Thickness-to-chord (\$\\tau\$)")
#     ax[end].set_xlabel("\$c_l\$")
#     ax[1].set_ylabel("\$c_d\$")

#     ax[2].plot(airf.cl, airf.A[end, :, :, 3], label = airf.τ)
#     ax[2].set_ylabel("\$c_m\$")
# end
