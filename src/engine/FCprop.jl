using PyPlot
include("gasfun.jl")
include("gascalc.jl")
include("hxfun.jl")
include("PEMfuelcell.jl")

#Parameters that can be optimized
η_p = 0.7
type = "HT-PEMFC"
#Set HX targets
ε = 0.5

#Design conditions: takeoff
P = 170 * 15e3 #takeoff stack power

Ta = 303.15 #K, corresponding to 30 degrees C
Mp_in  = 0.1
mdot_p = 40


HXgas_NaN = HX_gas("0","0", [NaN], NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
HXgeom_NaN = HX_tubular(0, 0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "")

HXgas = deepcopy(HXgas_NaN)
HXgeom = deepcopy(HXgeom_NaN)
#---------------------------------
# PEMFC calculations
#---------------------------------
p0 = 101325
V_stack = 300
# Inputs
p_A = 3 * p0
p_C = 3 * p0
λ_H2 = 3
λ_O2 = 3
t_M = 25e-6
t_A = 350e-6
t_C = 350e-6

u = PEMFC_inputs(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "")
u.p_A = p_A
u.p_C = p_C
u.λ_H2 = λ_H2
u.λ_O2 = λ_O2
u.t_M = t_M
u.t_A = t_A
u.t_C = t_C
u.type = type

#Find cell voltage depending on type

if type == "LT-PEMFC"
    u.j = 1.6e4
    u.T = 353.15
    u.x_H2O_A = water_sat_pressure(u.T) / p_A #fully saturated gas
    u.x_H2O_C = water_sat_pressure(u.T) / p_C #fully saturated gas

else
    u.j = 2.0e4
    u.T = 453.15
    u.x_H2O_A = 0.2#water_sat_pressure(T) / p_A #fully saturated gas
    u.x_H2O_C = 0.2#water_sat_pressure(T) / p_C #fully saturated gas

end

#PEMFC calculations
n_cells, A_cell, Q = PEMsize(P, V_stack, u)

#V_od, Q_od = PEMoper(Pi, n_cells, A_cell, u)

#---------------------------------
# Design HX for takeoff conditions
#---------------------------------
#TODO: ideally want to do multipoint design, where we design for different mission points and choose largest HX
#Fluid parameters
Tp_in = Ta
Tc_in = u.T
pp_in = p0
pc_in = p_A
fluid_p = "air"
alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]

if type == "LT-PEMFC"
    fluid_c = "liquid water"
else
    fluid_c = "liquid ethylene glycol"
end

#specific heats
_, cp_c, _, _, _, _ = liquid_properties(fluid_c, Tc_in)
_, _, _, _, cp_p, Rp = gassum(alpha_p, length(alpha_p), Tp_in)

#Calculate minimum heat capacity rate from desired effectiveness and temperature difference
C_min = abs(Q / (ε * (Tp_in - Tc_in)))

#Design for C_min being C_c
mdot_c = C_min / cp_c

ρ_p = pp_in / (Rp * Tp_in)
γ = cp_p / (cp_p - Rp)
A_cs = mdot_p / (ρ_p * Mp_in * sqrt(γ * Rp * Tp_in))
l = sqrt(A_cs)

HXgas.fluid_p = fluid_p
HXgas.fluid_c = fluid_c
HXgas.mdot_p = mdot_p
HXgas.mdot_c = mdot_c
HXgas.ε = ε
HXgas.Tp_in = Tp_in
HXgas.Tc_in = Tc_in
HXgas.Mp_in  = Mp_in
HXgas.pp_in = pp_in
HXgas.pc_in = pc_in

HXgas.alpha_p = alpha_p

HXgeom.fconc = 0
HXgeom.frecirc = 0
HXgeom.l = l
HXgeom.xl_D = 1
HXgeom.Rfp = 0.0001763
HXgeom.Rfc = 0.0001763
HXgeom.material = "SS304"

initial_x = [0.1, 6, 4]

hxoptim!(HXgas, HXgeom, initial_x)
hxsize!(HXgas, HXgeom)

#---------------------------------
# HX in off-design conditions
#---------------------------------

HXgas_od = deepcopy(HXgas_NaN)

#Off design conditions: takeoff at lower power
P_od = 0.7 * P #takeoff stack power

mdot_pod = mdot_p * P_od / P

_, Q_od = PEMoper(P_od, n_cells, A_cell, u)

HXgas_od.fluid_p = fluid_p
HXgas_od.fluid_c = fluid_c
HXgas_od.mdot_p = mdot_pod
HXgas_od.Tp_in = Tp_in
HXgas_od.Tc_in = Tc_in
HXgas_od.pp_in = pp_in
HXgas_od.pc_in = pc_in
HXgas_od.alpha_p = alpha_p

HXoffDesignCalc!(HXgas_od, HXgeom, Q_od)