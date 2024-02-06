#include("../src/structures/structures.jl")
include("../src/structures/tanksize.jl")
include("../src/structures/tankWmech.jl")
include("../src/structures/tankWthermal.jl")
using Roots
using NLsolve

function insulation_conductivity_calc(T, material)
    if material == "polyurethane"
        k = 5e-3 + 20e-3/180 * (T - 20) #Linear fit to Fig. 4.78 in Brewer (1991)
    end
    return k
end

mutable struct thermal_params
    l_cyl::Float64
    l_tank::Float64
    r_tank::Float64
    Shead::Array{Float64}
    hconvgas::Float64
    hconvair::Float64
    t_cond::Array{Float64} 
    material::String
    Tfuel::Float64
    Tair::Float64
    ifuel::Int64
    thermal_params() = new() 
end

function residuals_Q(x, p)

    Q = x[1]
    T_w = x[2]
    T_mli = x[3:end]

    #Unpack parameters
    l_cyl = p.l_cyl
    l_tank = p.l_tank
    r_tank = p.r_tank
    Shead = p.Shead
    hconvgas = p.hconvgas
    hconvair = p.hconvair
    t_cond = p.t_cond
    material = p.material
    Tfuel = p.Tfuel
    Tair = p.Tair
    ifuel = p.ifuel

    r_inner = r_tank #- thickness
    ΔT = Tair - Tfuel
    thickness = sum(t_cond)  # total thickness of insulation

    # Radiation
    σ = 5.6704e-8 #W/(m^2 K^4), Stefan-Boltzmann constant
    ε = 0.95    # white aircraft (Verstraete)

    hradair = σ * ε * ((Tair^2) + (Tfuel^2)) * (Tair + Tfuel) #Radiative heat transfer coefficient; Eq. (2.28) in https://ahtt.mit.edu/
    h_air = hconvair + hradair # Combines radiative and convective heat transfer at outer end
    Rair_conv_rad = 1 / (h_air * (2π * r_tank * l_cyl + 2*Shead[end]))  # thermal resistance of ambient air (incl. conv and rad)

    S_int = (2π * (r_inner - thickness) * l_cyl) + 2*Shead[1] #liquid side surface area
    h_liq, _ = tank_heat_coeffs(T_w, ifuel, Tfuel, l_tank) #Find liquid-side heat transfer coefficient
    R_liq = 1 / (h_liq * S_int) #Liquid-side thermal resistance

    N = length(t_cond)       # Number of layers in insulation
    R_mli      = zeros(Float64, N)  #size of MLI resistance array (Based on number of layers)
    R_mli_ends = zeros(Float64, N)
    R_mli_cyl  = zeros(Float64, N)

    T_prev = T_w
    for i in 1:N
        k = insulation_conductivity_calc((T_mli[i] + T_prev)/2, material)
        R_mli_cyl[i] = log((r_inner  + t_cond[i])/ (r_inner)) / (2π*l_cyl * k) #Resistance of each MLI layer; from integration of Fourier's law in cylindrical coordinates
        R_mli_ends[i] = t_cond[i] / (k * 2*Shead[i])
        # Parallel addition of resistance
        R_mli[i] = (R_mli_ends[i] * R_mli_cyl[i]/(R_mli_ends[i] + R_mli_cyl[i])) 
        
        # Update r_inner
        r_inner = r_inner + t_cond[i]  
        T_prev = T_mli[i]
    end

    R_mli_tot = sum(R_mli)  #Total thermal resistance of MLI
    Req = R_mli_tot + R_liq + Rair_conv_rad # Total equivalent resistance of thermal circuit

    F = zeros(length(x))
    F[1] = Q - ΔT / Req

    T_calc = Tfuel + R_liq * Q
    F[2] = T_w - T_calc

    for i = 1:length(T_mli)
        T_calc = T_calc + R_mli[i] * Q
        F[i + 2] = T_mli[i] - T_calc
        
    end
    return F
end

const gee = 9.81
l_cyl = 2
l_tank = 4
r_tank = 1 
Shead = [1, 1, 1, 1, 1]
hconvgas= 0 
hconvair = 15
t_cond = [0.05, 1.524e-5, 0.05, 1.524e-5, 1.57e-2] #m
k = ones(length(t_cond)) * 5e-3
material = "polyurethane"

Tfuel = 20
Tair = 220
time_flight = 5*3600
ifuel = 40

p = thermal_params()
p.l_cyl = l_cyl
p.l_tank = l_tank
p.r_tank = r_tank
p.Shead = Shead
p.hconvgas = hconvgas
p.hconvair = hconvair
p.t_cond = t_cond
p.material = material
p.Tfuel = Tfuel
p.Tair = Tair
p.ifuel = ifuel

thickness = sum(t_cond)  # total thickness of insulation
ΔT = Tair - Tfuel
qfac = 1.3         # Account for heat leak from pipes and valves

fun(x) = residuals_Q(x, p)

guess = zeros(length(t_cond) + 2)
guess[1] = 100
guess[2] = Tfuel + 1

for i = 1:length(t_cond)
    guess[i + 2] = Tfuel + ΔT * sum(t_cond[1:i])/ thickness
end
sol = nlsolve(fun, guess)

_, h_v = tank_heat_coeffs(Tfuel, ifuel, Tfuel, l_tank) #Liquid side h and heat of vaporization

Q = qfac * sol.zero[1]    # Heat flux from ambient to cryo fuel, including extra heat leak from valves etc as in eq 3.20 by Verstraete
mdot_boiloff = Q / h_v  # Boil-off rate equals the heat flux divided by heat of vaporization
m_boiloff = mdot_boiloff * time_flight # Boil-off mass calculation





