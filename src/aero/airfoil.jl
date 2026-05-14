"""
    airfoil

A type representing a database of pre-computed airfoil data for a single Reynolds
number and a range of Mach numbers, sectional lift coefficients, and thickness-to-chord
ratios. By default, this is the C-series transonic airfoil family designed by M. Drela in MSES,
stored in `src/airfoil_data/C_airfoil.csv`. See [`airtable`](@ref) for the file schema and how the
struct is populated from CSV.

Overloads `Base.show` to print a summary of the `airfoil` model.

# Conventions and assumptions
- **4D data layout.** `A` is sized `(nMa, ncl, ntoc, nAfun)` where the 4th axis is **fixed-order**:
  `[CD, CDp, CDv, CDw, CM, alpha]` (i.e. `A[:,:,:,1]` is total `cd`, `A[:,:,:,6]` is `α`).
  This order is relied on by [`airfun`](@ref).
- **Single Reynolds number.** All data correspond to one fixed `Re`; chord-Re effects are
  applied externally as a power-law scaling (see [`wing_profiledrag_direct`](@ref)).
- **Units.** `Ma`, `cl`, `toc`, and aerodynamic coefficients are dimensionless. **`alpha`
  is stored in radians** (the source CSV is in degrees; conversion happens in [`airtable`](@ref)).
- **Infeasible (NaN) entries.** MSES-infeasible operating points (e.g. stalled, beyond drag
  divergence) are stored as `NaN` rather than duplicating neighbouring values (like original TASOPT). They are propagated
  by [`TASOPT.aerodynamics.spline`](@ref) and handled in [`airfun`](@ref) via CL-segment fallback + quadratic drag
  penalty. See [`airfoil_cl_limits`](@ref) to query the valid-data envelope.

# Fields
- `Ma::AbstractVector{Float64}`  : Mach numbers covered by the database.
- `cl::AbstractVector{Float64}`  : Sectional lift coefficients ".
- `toc::AbstractVector{Float64}` : Thickness-to-chord ratios ".
- `Re::Float64`                  : Reynolds number for which the data were computed.
- `A::AbstractArray{Float64}`    : Aerodynamic data, again, layout is [CD, CDp, CDv, CDw, CM, alpha].

Arrays of derivative arrays at the grid knots, used for tricubic interpolation in `airfun`.
- `A_M`        : `∂A/∂Ma`
- `A_toc`      : `∂A/∂toc`
- `A_cl`       : `∂A/∂cl`
- `A_M_toc`    : `∂²A/∂Ma∂toc`
- `A_M_cl`     : `∂²A/∂Ma∂cl`
- `A_cl_toc`   : `∂²A/∂cl∂toc`
- `A_M_cl_toc` : `∂³A/∂Ma∂cl∂toc`

See also [`airfun`](@ref), [`airtable`](@ref), [`airfoil_cl_limits`](@ref).
"""
struct airfoil{T<:AbstractFloat, 
        V<:AbstractVector{Float64}, 
        Ar<:AbstractArray{Float64}} 
    Ma::V
    cl::V
    toc::V
    Re::T # Data assumed for a single Re
  
    A::Ar # Airfoil aero data 
    
    A_M::Ar
    A_toc::Ar
    A_cl::Ar
    A_M_toc::Ar
    A_M_cl::Ar
    A_cl_toc::Ar
    A_M_cl_toc::Ar
end 

function Base.show(io::IO, airf::airfoil)
    printstyled(io, "Airfoil section database with:\n"; color=:bold)
    println(io, 
    "Re = ",airf.Re,"\n",
    "Ma = (", airf.Ma[1],", ", airf.Ma[end],")\n",
    "cl = (", airf.cl[1], ", ",airf.cl[end],")\n",
    "toc  = (", airf.toc[1], ", ",airf.toc[end],")\n"
    )
end

using Plots
"""
    plot_airf(airf::airfoil, iMach::Int = -1; use_interp::Bool=false, 
                    cl_range=range(minimum(airf.cl), maximum(airf.cl), length=100), 
                    toc_vals=airf.toc)

Plots aerodynamic characteristics of the airfoil database, including its drag and pitching moment curves across lift coefficients and airfoil thickness-to-chord ratios (`toc` = \$\\frac{t}{c}\$).  

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `airf::airfoil`: `airfoil` data struct 
    - `iMach::Int`: Index of Mach number to plot (default: -1, which selects the highest Mach number in the database).

    **Outputs:**
    - Returns a `Plots.Plot` object with three panels:
        - **Top subplot:** `cl` vs `cd` (drag polar for multiple `toc`).
        - **Bottom subplot:** `cl` vs `cm` (moment curve for multiple `toc`).
        - **Side panel:** Legend indicating the thickness-to-chord ratio values (`toc`).

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
                    toc_vals=airf.toc)
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
        n_toc = length(toc_vals)
        cd_total = zeros(n_cl, n_toc)
        cm_vals = zeros(n_cl, n_toc)
        aoa_vals = zeros(n_cl, n_toc)

        # Evaluate airfun for each combination of cl and toc
        for (j, toc) in enumerate(toc_vals)
            for (i, cl) in enumerate(cl_range)
                cdf, cdp, cdw, cm, aoa = airfun(cl, toc, mach, airf)
                cd_total[i, j] = cdf + cdp # Total drag (friction + pressure)
                cm_vals[i, j] = cm
                aoa_vals[i, j] = aoa
            end
        end
        cls = cl_range
        cds = cd_total
        cms = cm_vals
        aoas = aoa_vals
    else #just show the data points
        cls = airf.cl
        cds = airf.A[iMach,:,:,1]
        cms = airf.A[iMach,:,:,5]
        aoas = airf.A[iMach,:,:,6]
    end

    # Create two subplots
    p1 = plot(
        cls,
        cds,
        label = string.(airf.toc'),
        # xlabel = "\$c_l\$",
        ylabel = "\$c_d\$",
        legend=:false,
        grid = true,
    )
    
    p2 = plot(
        cls,
        cms,
        label = "toc = ".*string.(airf.toc'),
        ylabel = "\$c_m\$",
        grid = true,
        legend=:false,
    )
    
    p3 = plot(
        cls,
        rad2deg.(aoas),
        label = "toc = ".*string.(airf.toc'),
        xlabel = "\$c_l\$",
        ylabel = "\$α_⟂, \\mathrm{aoa}_⟂\$ [°]",
        grid = true,    
        legend=:false,
    )
    
    # If interp, overlay original data points as scatter
    if use_interp
        cl_scatter = airf.cl
        for (j, toc) in enumerate(airf.toc)
            cd_scatter = airf.A[iMach, :, j, 1]
            cm_scatter = airf.A[iMach, :, j, 5]
            aoa_scatter = rad2deg.(airf.A[iMach, :, j, 6])
            scatter!(p1, cl_scatter, cd_scatter; label="", marker=:circle, ms=3, color=j)
            scatter!(p2, cl_scatter, cm_scatter; label="", marker=:circle, ms=3, color=j)
            scatter!(p3, cl_scatter, aoa_scatter; label="", marker=:circle, ms=3, color=j)
        end
    end

    #legend
    labels = string.(airf.toc')
    p_legend = plot((1:length(airf.toc))', labels = labels,
        legendtitle = "Thickness-to-chord \n(toc)",
        legend_title_font_pointsize = 7,
        legendfontsize=7, legend=:outertop, legendcolumns=1,
        fg_color_legend = nothing, frame=:none)
    
    l = @layout [[a; b; c] d{0.2w}]
    # Combine the subplots vertically
    f1 = plot(p1, p2, p3, p_legend, layout = l, link = :x,
        size=(600, 600),
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
