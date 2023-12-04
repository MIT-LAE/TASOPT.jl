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
HXgas_NaN = HX_gas("0","0", [NaN], NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
HXgeom_NaN = HX_tubular(0, 0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)

#---------------------------------     
# hxsize!()
#---------------------------------
HXgas = deepcopy(HXgas_NaN)
HXgeom = deepcopy(HXgeom_NaN)

HXgas.gas_h = "air"
HXgas.gas_c = "h2"
HXgas.mdot_h = 1144/60
HXgas.mdot_c = 9.95/60
HXgas.ε = 0.8
HXgas.Th_in = 778
HXgas.Tc_in = 264
HXgas.Mh_in  = 0.19
HXgas.Mc_in = 0.0285
HXgas.ph_in = 40e3
HXgas.pc_in = 1515e3

HXgas.alpha_h = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
HXgas.igas_c = 40

HXgeom.fconc = 1
HXgeom.frecirc = 0
HXgeom.D_i = 0.564
HXgeom.t = 0.03e-2 #m, wall thicknesss
HXgeom.l = 0.6084530646014857 #tube length
HXgeom.n_stages = 4
HXgeom.xt_D = 6
HXgeom.xl_D = 1.25
HXgeom.kw = 45 #thermal conductivity of steel, W/m/K
HXgeom.Rfh = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
HXgeom.ρw = 7930

hxsize!(HXgas, HXgeom)

size_out = [HXgas.Th_out, HXgas.Tc_out, HXgas.Δp_h, HXgeom.N_t, HXgeom.n_passes, HXgeom.tD_o, HXgeom.A_cs]

size_out_check = 
[ 731.5893294226444, 665.8919656336637, 1432.24788197369, 62.03560518812884,  7.9999971021522125,
  0.004760326082769499, 1.0189779296746375]

size_check = maximum(abs.(size_out - size_out_check)./size_out_check) #Find maximum relative difference

if size_check < 1e-3
    println("hxsize!() passes test")
else
    println("hxsize!() FAILS test")
end

#---------------------------------     
# hxweight()
#---------------------------------
gee = 9.81
fouter = 1

W = hxweight(gee, HXgeom, fouter)

W_check = 790.0013892027117

if abs(W-W_check) / W_check < 1e-4
    println("hxweight() passes test")
else
    println("hxweight() FAILS test")
end

#---------------------------------     
# hxoper!()
#---------------------------------
HXgas = deepcopy(HXgas_NaN)
HXgeom = deepcopy(HXgeom_NaN)

HXgas.gas_h = "air"
HXgas.gas_c = "h2"
HXgas.mdot_h = 2*1144/60
HXgas.mdot_c = 2*9.95/60
HXgas.ε = 0.8
HXgas.Th_in = 778
HXgas.Tc_in = 264
HXgas.Mh_in  = 0.19
HXgas.Mc_in = 0.0285
HXgas.ph_in = 3*40e3
HXgas.pc_in = 3*1515e3

HXgas.alpha_h = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
HXgas.igas_c = 40

HXgeom.fconc = 1
HXgeom.frecirc = 0
HXgeom.D_i = 0.564
HXgeom.t = 0.03e-2 #m, wall thicknesss
HXgeom.tD_o = 0.004760326082769499
HXgeom.A_cs = 1.0189779296746375
HXgeom.l = 0.6084530646014857 #tube length
HXgeom.n_stages = 4
HXgeom.n_passes = 8
HXgeom.N_t = 62
HXgeom.xt_D = 6
HXgeom.xl_D = 1.25
HXgeom.kw = 45 #thermal conductivity of steel, W/m/K
HXgeom.Rfh = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
HXgeom.ρw = 7930

hxoper!(HXgas, HXgeom)

oper_out = [HXgas.Th_out, HXgas.Tc_out, HXgas.Δp_h, HXgas.ε]

oper_out_check = [740.2023159693746, 592.0560174320444, 1735.5866896218854, 0.6521339168630342]

oper_check = maximum(abs.(oper_out - oper_out_check)./oper_out_check) #Find maximum relative difference

if oper_check < 1e-3
    println("hxoper!() passes test")
else
    println("hxoper!() FAILS test")
end

#---------------------------------     
# hxoptim!()
#---------------------------------
HXgas = deepcopy(HXgas_NaN)
HXgeom = deepcopy(HXgeom_NaN)

HXgas.gas_h = "air"
HXgas.gas_c = "h2"
HXgas.mdot_h = 1144/60
HXgas.mdot_c = 9.95/60
HXgas.ε = 0.8
HXgas.Th_in = 778
HXgas.Tc_in = 264
HXgas.Mh_in  = 0.19
HXgas.ph_in = 40e3
HXgas.pc_in = 1515e3

HXgas.alpha_h = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
HXgas.igas_c = 40

HXgeom.fconc = 1
HXgeom.frecirc = 0
HXgeom.D_i = 0.564
HXgeom.t = 0.03e-2 #m, wall thicknesss
HXgeom.l = 0.6084530646014857 #tube length
HXgeom.xl_D = 1
HXgeom.kw = 45 #thermal conductivity of steel, W/m/K
HXgeom.Rfh = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
HXgeom.ρw = 7930

#Calculate starting point
#First calculate minimum tube length
_, _, _, _, cp_h_in, Rh = gassum(HXgas.alpha_h, length(HXgas.alpha_h), HXgas.Th_in)
γ_h_in = cp_h_in / (cp_h_in - Rh)
ρ_h_in = HXgas.ph_in / (Rh * HXgas.Th_in)
Vh_in = HXgas.Mh_in * sqrt(γ_h_in * Rh * HXgas.Th_in)

A_cs = HXgas.mdot_h / (ρ_h_in * Vh_in) #Cross-sectional area of freestream

if HXgeom.fconc #Flow is concentric
    D_i = HXgeom.D_i
    D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

    lmin = (D_o - D_i) / 2 #minimum tube length
    linit = 1.1 * lmin
else #square cross-section
    AR_min = 0.1 #Minimum aspect ratio
    lmin = sqrt(AR_min * A_cs)
    linit = sqrt(A_cs)
end

#Now set starting point
initial_x = [3, linit, 4, 4] #Initial guess
 
# Optimize heat exchanger design parameters
hxoptim!(HXgas, HXgeom, initial_x)

optim_out = [HXgas.Mc_in, HXgeom.l, HXgeom.n_stages, HXgeom.xt_D]

hxsize!(HXgas, HXgeom)

I = HXgas.Pl_h + HXgas.Pl_c #Optimizer may choose slightly different points with similar objective function. Check I instead of optim outputs

I_check = 73945.57024509525

if abs(I-I_check) / I_check < 1e-2
    println("hxoptim!() passes test")
else
    println("hxoptim!() FAILS test")
end

#---------------------------------     
# hxoptim!() with recirculation and rectangular
#---------------------------------
HXgas = deepcopy(HXgas_NaN)
HXgeom = deepcopy(HXgeom_NaN)

HXgas.gas_h = "air"
HXgas.gas_c = "h2"
HXgas.mdot_h = 2 * 49.9/60
HXgas.mdot_c = 9.95/60
HXgas.ε = 0.825
HXgas.Th_in = 791
HXgas.Tc_in = 20
HXgas.Mh_in  = 0.01
HXgas.ph_in = 1515e3
HXgas.pc_in = 1515e3
HXgas.recircT  = 200
HXgas.h_lat  = 446e3 + 670e3

HXgas.alpha_h = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
HXgas.igas_c = 40

HXgeom.fconc = 0
HXgeom.frecirc = 1
HXgeom.t = 0.03e-2 #m, wall thicknesss
HXgeom.xl_D = 1
HXgeom.kw = 45 #thermal conductivity of steel, W/m/K
HXgeom.Rfh = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
HXgeom.ρw = 7930

#Calculate starting point
#First calculate minimum tube length
_, _, _, _, cp_h_in, Rh = gassum(HXgas.alpha_h, length(HXgas.alpha_h), HXgas.Th_in)
γ_h_in = cp_h_in / (cp_h_in - Rh)
ρ_h_in = HXgas.ph_in / (Rh * HXgas.Th_in)
Vh_in = HXgas.Mh_in * sqrt(γ_h_in * Rh * HXgas.Th_in)

A_cs = HXgas.mdot_h / (ρ_h_in * Vh_in) #Cross-sectional area of freestream

if HXgeom.fconc #Flow is concentric
    D_i = HXgeom.D_i
    D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

    lmin = (D_o - D_i) / 2 #minimum tube length
    linit = 1.1 * lmin
else #square cross-section
    AR_min = 0.1 #Minimum aspect ratio
    lmin = sqrt(AR_min * A_cs)
    linit = sqrt(A_cs)
end

#Now set starting point
initial_x = [3, linit, 4, 4] #Initial guess

hxoptim!(HXgas, HXgeom, initial_x)

hxsize!(HXgas, HXgeom)

I = HXgas.Pl_h + HXgas.Pl_c #Optimizer may choose slightly different points with similar objective function. Check I instead of optim outputs

optimrec_out = [HXgas.Mc_in, HXgas.mdot_r, HXgeom.l, HXgeom.n_stages, HXgeom.xt_D]

I = HXgas.Pl_h + HXgas.Pl_c #Optimizer may choose slightly different points with similar objective function. Check I instead of optim outputs

I_check = 2681.660392262237

if abs(I-I_check) / I_check < 1e-2
    println("hxoptim!() with recirculation passes test")
else
    println("hxoptim!() with recirculation FAILS test")
end

# println(Mc_in,'|', n_passes,'|', n_stages, '|',xt_D,'|', xl_D)
# println(N_t,'|',l,"|",Δp_h/ph_in)