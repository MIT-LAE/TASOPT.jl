#----------------------------------------------------------------
#     Heat exchanger sizing tests
#
#     Thos script determines the geometric parameters for a cross-flow heat 
#     exchanger in air. It assumes that air in flowing through the exchanger
#     and heat is absorbed by a coolant fluid.

using BenchmarkTools
#---------------------------------     
# Preamble
#---------------------------------
include("../src/engine/hxfun.jl")
include("../src/engine/gascalc.jl")
include("../src/engine/gasfun.jl")

#---------------------------------     
# Test functions individually
#---------------------------------

#---------------------------------     
# hxsize()
#---------------------------------
gas_h = "air"
gas_c = "h2"
mdot_h = 1144/60
mdot_c = 9.95/60
ε = 0.8
Th_in = 778
Tc_in = 264
Mh_in  = 0.19
Mc_in = 0.0285
ph_in = 40e3
pc_in = 1515e3
D_i = 0.564
t = 0.03e-2 #m, wall thicknesss
n_passes = 8
n_stages = 4
xt_D = 6
xl_D = 1.25
kw = 45 #thermal conductivity of steel, W/m/K
Rfh = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W

size_out = hxsize(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, Mc_in, ph_in, pc_in, D_i, t, 
n_passes, n_stages, xt_D, xl_D, kw, Rfh, Rfc)

size_out = [i for i in size_out]

size_out_check = 
[731.6442346709227, 675.2, 1446.9828110715948, 0.004763196129933139, 
61.99822585125962, 0.611259573256122, 1.0189779296746375, 1.271022814615896, 149669.71582976307, 10752.604499917124]

size_check = maximum(abs.(size_out - size_out_check)./size_out_check) #Find maximum relative difference

if size_check < 1e-3
    println("hxsize() passes test")
else
    println("hxsize() FAILS test")
end

#---------------------------------     
# hxoper()
#---------------------------------
gas_h = "air"
gas_c = "h2"
mdot_h = 2*1144/60
mdot_c = 2*9.95/60
Th_in = 778
Tc_in = 264
ph_in = 3*40e3
pc_in = 3*1515e3
A_cs = 1.0189779296746375
tD_o = 0.004763196129933139
t = 0.03e-2 #m, wall thicknesss
l = 0.611259573256122
N_t = 61.99822585125962
n_passes = 8
n_stages = 4
xl_D = 1.25
kw = 45 #thermal conductivity of steel, W/m/K
Rfh = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W

oper_out = hxoper(gas_h, gas_c, mdot_h, mdot_c, Th_in, Tc_in, ph_in, pc_in, A_cs, tD_o, t, 
l, N_t, n_passes, n_stages, xl_D, kw, Rfh, Rfc)

oper_out = [i for i in oper_out]

oper_out_check = [740.1945947041926, 599.3538129999303, 1755.6285075738906]

oper_check = maximum(abs.(oper_out - oper_out_check)./oper_out_check) #Find maximum relative difference

if oper_check < 1e-3
    println("hxoper() passes test")
else
    println("hxoper() FAILS test")
end

#---------------------------------     
# hxoptim()
#---------------------------------
gas_h = "air"
gas_c = "h2"
mdot_h = 1144/60
mdot_c = 9.95/60
ε = 0.8
Th_in = 778
Tc_in = 264
Mh_in  = 0.19
ph_in = 40e3
pc_in = 1515e3
D_i = 0.564
t = 0.03e-2 #m, wall thicknesss
xl_D = 1
kw = 45 #thermal conductivity of steel, W/m/K
Rfh = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W

# Optimize heat exchanger design parameters
Mc_in, n_passes, n_stages, xt_D = hxoptim(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, ph_in, pc_in, 
D_i, t, xl_D, kw, Rfh, Rfc)

optim_out = [Mc_in, n_passes, n_stages, xt_D]

_, _, _, _, _, _, _, _, Pl_h, Pl_c = hxsize(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, Mc_in, ph_in, pc_in, D_i, t, 
n_passes, n_stages, xt_D, xl_D, kw, Rfh, Rfc)

I = Pl_h + Pl_c #Optimizer may choose slightly different points with similar objective function. Check I instead of optim outputs

I_check = 75782.45740757471

if abs(I-I_check) / I_check < 1e-2
    println("hxoptim() passes test")
else
    println("hxoptim() FAILS test")
end

#---------------------------------     
# hxweight()
#---------------------------------
gee = 9.81
fouter = 1
ρ_m = 7930 #density of steel (kg/m^3)
tD_o = 0.004763196129933139
t = 0.03e-2 #m, wall thicknesss
l = 0.611259573256122
N_t = 61.99822585125962
n_passes = 8
n_stages = 4

W = hxweight(gee, fouter, ρ_m, tD_o, t, l, N_t, n_passes, n_stages)

W_check = 793.6777467557885

if abs(W-W_check) / W_check < 1e-4
    println("hxweight() passes test")
else
    println("hxweight() FAILS test")
end

# println(Mc_in,'|', n_passes,'|', n_stages, '|',xt_D,'|', xl_D)
# println(N_t,'|',l,"|",Δp_h/ph_in)