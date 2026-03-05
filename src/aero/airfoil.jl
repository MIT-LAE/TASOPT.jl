"""
    airfoil

A type representing a database of pre-computed airfoil data for a single Reyolds number 
and a range of Mach numbers, sectional lift coefficients, and thickness-to-chord ratios.
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
struct airfoil{T<:AbstractFloat, 
        V<:AbstractVector{Float64}, 
        Ar<:AbstractArray{Float64}} 
    Ma::V
    cl::V
    τ::V
    Re::T # Data assumed for a single Re
  
    A::Ar # Airfoil aero data 
    
    A_M::Ar
    A_τ::Ar
    A_cl::Ar
    A_M_τ::Ar
    A_M_cl::Ar
    A_cl_τ::Ar
    A_M_cl_τ::Ar
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

using Plots
"""
    plot_airf(airf::airfoil, iMach::Int=-1)

Plots aerodynamic characteristics of the airfoil database, including its drag and pitching moment curves across lift coefficients and airfoil thickness-to-chord ratios (`τ`).  

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `airf::airfoil`: `airfoil` data struct 
    - `iMach::Int`: Index of Mach number to plot (default: -1, which selects the highest Mach number in the database).

    **Outputs:**
    - Returns a `Plots.Plot` object with three panels:
        - **Top subplot:** `cl` vs `cd` (drag polar for multiple `τ`).
        - **Bottom subplot:** `cl` vs `cm` (moment curve for multiple `τ`).
        - **Side panel:** Legend indicating the thickness-to-chord ratio values (`τ`).

    Note: can also be called with `iMach::Colon` to generate plots for all Mach numbers in the database, or with an aircraft object to plot its wing airfoil data directly.

    Sample usage:

        ```julia
        ac = load_default_model(); size_aircraft!(ac)
        airf = ac.wing.airsection
        f1 = plot_airf(airf, 1)
        f2 = plot_airf(ac)
        ```
"""
function plot_airf(airf::airfoil, iMach::Int = -1; use_interp::Bool=false, 
                    cl_range=range(minimum(airf.cl), maximum(airf.cl), length=100), 
                    tau_vals=airf.τ)
    #mach number for plotting
    if iMach == -1 
        iMach = length(airf.Ma) 
    end
    mach = airf.Ma[iMach]

    #get data
    #if plotting the tricubic interpolation
    if use_interp
        # Pre-allocate arrays for refined data
        n_cl = length(cl_range)
        n_tau = length(tau_vals)
        cd_total = zeros(n_cl, n_tau)
        cm_vals = zeros(n_cl, n_tau)

        # Evaluate airfun for each combination of cl and tau
        for (j, tau) in enumerate(tau_vals)
            for (i, cl) in enumerate(cl_range)
                cdf, cdp, cdw, cm, alpha = airfun(cl, tau, mach, airf)
                cd_total[i, j] = cdf + cdp # Total drag (friction + pressure)
                cm_vals[i, j] = cm
            end
        end
        cls = cl_range
        cds = cd_total
        cms = cm_vals
    else #just show the data points
        cls = airf.cl
        cds = airf.A[iMach,:,:,1]
        cms = airf.A[iMach,:,:,5]
    end

    # Create two subplots
    p1 = plot(
        cls,
        cds,
        label = string.(airf.τ'),
        # xlabel = "\$c_l\$",
        ylabel = "\$c_d\$",
        legend=:false,
        grid = true,
    )
    
    p2 = plot(
        cls,
        cms,
        label = "τ = ".*string.(airf.τ'),
        xlabel = "\$c_l\$",
        ylabel = "\$c_m\$",
        grid = true,
        legend=:false,
    )
    
    # If interp, overlay original data points as scatter
    if use_interp
        cl_scatter = airf.cl
        for (j, tau) in enumerate(airf.τ)
            cd_scatter = airf.A[iMach, :, j, 1]
            cm_scatter = airf.A[iMach, :, j, 5]
            scatter!(p1, cl_scatter, cd_scatter; label="", marker=:circle, ms=3, color=j)
            scatter!(p2, cl_scatter, cm_scatter; label="", marker=:circle, ms=3, color=j)
        end
    end

    #legend
    labels = string.(airf.τ')
    p_legend = plot((1:length(airf.τ))', labels = labels,
        legendtitle = "Thickness-to-chord (τ)",
        legend_title_font_pointsize = 7,
        legendfontsize=7, legend=:outertop, legendcolumns=1,
        fg_color_legend = nothing, frame=:none)
    
    l = @layout [[a; b] c{0.2w}]
    # Combine the subplots vertically
    f1 = plot(p1, p2, p_legend, layout = l, link = :x,
        suptitle="Airfoil Section Database, Mach = $(mach)")
    return f1
end

function plot_airf(airf::airfoil, iMach::Colon; use_interp::Bool=false)
    nMa = length(airf.Ma)
    figs = []
    for iMa in 1:nMa
        fig = plot_airf(airf, iMa,use_interp=use_interp)
        push!(figs, fig)
        display(fig)
    end
    return figs
end
