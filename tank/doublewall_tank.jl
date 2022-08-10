using Base.Iterators
using Roots
using Printf
using PyCall
pygui()
using PyPlot
pygui(true)


"""
doublewalled_tank calcualtes the weight and heat transfer of a 
double walled LH2 tank

## Inputs
fuse_clearance: [m] Clearance from fuselage
d_vaccuum     : [m] distance between inner and outer vessel
Ltank         : [m] Outer most length of tank

p             : [Pa] Internal storage pressure
σa_inner      : [Pa] Allowable stress in inner vessel ~129.2 MPa for SS 316
ew            : [-]  Weld efficiency - Fully radiographed -> 1.0 (See Barron pg 361)
AR            : [-]  Aspect ratio of heads - 1 => hemisphere
ρinner        : [kg/m³] Density of inner tank material. High Nickel steel ideal to avoid embrittilement (SS314 = 7900 kg/m³)
Nstiff_in     : [-] Number of inner vessel stiffners
θin_support   : [rad] Location of inner vessel support 

σa_outer      : [Pa] Allowable stress in outer vessel
poiss         : [-] Poission ratio of outer vessel material
E_outer       : [Pa] Young's modulus of outer vessel material
ρouter        : [kg/m³] Density of outer tank material. 

ullage        : [-]   Ullage fraction - empty space in the tank for vapor
ρfuel         : [kg/m³] fuel density

"""
function doublewalled_tank(Rfuse::Float64, dRfuse::Float64, fuse_clearance, d_vaccuum,  Ltank,
    p::Float64, σa_inner::Float64, ew::Float64, AR, ρinner, Nstiff_in, θin_support,
     σa_outer, poiss, E_outer, ρouter,
     ullage, ρfuel,)

    Roo = Rfuse - fuse_clearance  # Outer Radius of outer vessel
    Doo = 2*Roo                   # OD of outer vessel

    Roi = Roo - d_vaccuum         # Outer Radius of inner vessel

    L = Ltank - 2*d_vaccuum - 2*Roo/AR # Length of cylindrical portion of inner tank

    Lh_inner = Roi/AR
    Lh_outer = Lh_inner + d_vaccuum # Assume inner and outer vessels have same cylindrical length

    # Inner vessel
        Doi = 2*Roi # OD of inner vessel
        t_cyl  = p*Doi/(2*σa_inner*ew + 0.8*p) # Min thickness of cylindrical vessel ASME pressure vessel
        # println("t_cyl = $t_cyl")
        K = (1/6) * (AR^2 + 2)            # Aspect ratio of 2:1 for the head (# Barron pg 359)
        t_head = p*Doi*K/(2*σa_inner*ew + 2*p*(K-1)) # Min thickness of hemisphirical/ ellpitical heads ASME pressure vessel
        # println("t_head = $t_head")
        
        ## Areas
        Ahead = 2π*Roi^2 * (1.0/3.0 + 2.0/3.0*(1/AR)^1.6) ^ (1/1.6) # Surface area of head
        Acyl  = 2π*Roi*t_cyl  # Cross-sectional area
        
        ## Volume and Weight
        Vcyl  = Acyl*L
        Vhead = Ahead*t_head
        
        Wcyl  = Vcyl*ρinner * gee
        Whead = Vhead*ρinner * gee
        Wtank_inner = Wcyl + Whead
        println("Wtank_inner = $Wtank_inner")
        
        ## Fuel Volume
        Rii  = Roi - t_cyl # internal radius
        Rih = Roi - t_head
        V_ellipsoid = 2π*(Rih^3/AR)/3
        V_cylinder  = π*Rii^2*L
        V_internal  = V_cylinder + V_ellipsoid
        V_fuel = V_internal/(1 + ullage)

        Wfuel_tot = V_fuel * ρfuel * gee
        println("Wfuel_tot = $Wfuel_tot")
        W_inner = Wfuel_tot + Wtank_inner
        
        ## Stiffening
        ϕlist, klist, ϕmax, kmax = stiffeners_bendingM(θin_support)
        Mmax = kmax * W_inner/Nstiff_in * Roi/(2π) #use outer dia for now
        Zmax = Mmax/σa_inner
        
        #Use standard metric beams to get weights per unit length
        if Zmax < 89.9e-6
            # W 100 x 100 x 19.3
            W′ = 19.3 # kg/m
        elseif Zmax < 139.5e-6
            # W 130 x 130 x 23.8
            W′ = 23.8 # kg/m
        end
        Perimeter = 2π*Roi
        W_stiffners = Nstiff_in * Perimeter * W′ * gee
        println("W_stiffners = $W_stiffners")
        # Outer vessel
        # Collapsing pressure for "long" cylinder
        # pc = 2*E(t_cyl/Doo)^3 / (1 - poiss^2)
        pc = 4*pSL
        t_cyl = Doo*cbrt(pc*(1-poiss^2)/(2*E_outer))
        # println("t_cyl = $t_cyl")
        
        if L/Doo < 1.140*(1-poiss^2)^(1.0/4)*sqrt(Doo/t_cyl) 
            # println("Warning: Outer Cycl isn't really 'long' - using 'short' tank relation instead")
            t_cyl = tank_buckling(pc, E_outer, poiss, L, Doo, t_cyl)
            # println("t_cyl = $t_cyl")
            L/Doo > 1.140*(1-poiss^2)^(1.0/4)*sqrt(Doo/t_cyl) && println("ohhh snap")
        end
        
        K1 = 0.90 # See table 7.6 in Barron p. 367
        t_head = K1*Doo*sqrt(pc*sqrt(3*(1 - poiss^2))/ (0.5*E_outer))
        # println("t_head = $t_head")
        
        ## Areas
        Ahead = 2π*Roo^2 * (1.0/3.0 + 2.0/3.0*(Lh_outer/Roo)^1.6) ^ (1/1.6) # Surface area of head
        Acyl  = 2π*Roo*t_cyl  # Cross-sectional area

        ## Volume and Weight
        Vcyl  = Acyl*L
        Vhead = Ahead*t_head

        Wcyl  = Vcyl *ρouter*gee
        Whead = Vhead*ρouter*gee   
        Wtank_outer = Wcyl + Whead
        println("Wtank_outer = $Wtank_outer")

        tank_total = Wtank_inner + Wtank_outer + W_stiffners + Wfuel_tot
        η = Wfuel_tot/ tank_total *100
        println("Tank total = $tank_total, η = $η")
        println("fadd = ", W_stiffners/(tank_total - W_stiffners))
        println((tank_total - W_stiffners)*(1+0.032413))

end

"""
## Stiffening rings
returns 2πM/(WR)
"""
function stiffeners_bendingM(θ)
    ϕlist = LinRange(0.0, π, 1800)
    k = zeros(length(ϕlist))

    for (i,ϕ) in enumerate(ϕlist)
        if 0 ≤ ϕ ≤ θ
            k[i] = 0.5*cos(ϕ) + ϕ*sin(ϕ) - (π - θ)*sin(θ) + cos(θ) + cos(ϕ)*(sin(θ)^2)
        elseif θ ≤ ϕ ≤ π
            k[i] = 0.5*cos(ϕ) - (π - ϕ)*sin(ϕ) + θ*sin(θ)  + cos(θ) + cos(ϕ)*(sin(θ)^2)
        end
    end
    kmax, imax = findmax(abs.(k))
    ϕmax = ϕlist[imax]
    return ϕlist, k, ϕmax, kmax
end

"""
Visualize the bending moment distribution
"""
function plot_stiffeners()
    plt.style.use(["./prash.mplstyle", "tableau-colorblind10"])
    fig, ax = plt.subplots(figsize = (8,5))
    for θ in [30, 60, 80]
        ϕ, k, ϕmax, kmax = stiffeners_bendingM(θ*π/180.0)
        
        ax.plot(ϕ.*180.0./π, k, label = "\$\\theta \$ = $θ\$^\\circ, \\phi_{max} = \$ $(round(ϕmax*180.0/π; digits = 1))")
        ax.axvline(ϕmax*180.0/π, lw = 1.0, alpha = 0.5, color = "k")
    end
    ax.legend()
    ax.set_xlabel("Location \$\\phi\$")
    ax.set_ylabel("\$\\frac{2\\pi M}{WR}\$")
end

function tank_buckling(pc, E, poiss, L, Do, tguess)
    f(t_D) = 2.42*E*(t_D)^(5.0/2) / ( (1-poiss^2)^(3.0/4.0) * ((L/Do) - 0.45*sqrt(t_D)) ) - pc
    t_D = find_zero(f, tguess)
    t = t_D*Do
end

gee = 9.81
pSL = 101325.0

Rfuse = 2.9
dRfuse = 0.0
fuse_clearance = 0.05
d_vaccuum = 0.1
Ltank = 7.0
p = 2*101.325e3
σa_inner = 130.0e6
ew = 1.0
AR = 2.0
ρinner = 2700.0#7900.0
Nstiff_in = 3
θin_support = 80*π/180.0
σa_outer = 130.0e6
poiss = 0.33
E_outer = 69.0e9
ρouter = 2700.0
ullage = 0.1
ρfuel = 67.0
doublewalled_tank(Rfuse, dRfuse, fuse_clearance, d_vaccuum,  Ltank,
    p, σa_inner, ew, AR, ρinner, Nstiff_in, θin_support,
     σa_outer, poiss, E_outer, ρouter,
     ullage, ρfuel,)