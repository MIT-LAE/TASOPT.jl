using TASOPT
using NLsolve
using Roots

mutable struct innertank()
    ptank::Float64
    sigskin::Float64
    rhoskintank::Float64
    ARtank::Float64
    clearance_fuse::Float64
    material_insul::Vector{String}
    ftankstiff::Float64
    ftankadd::Float64
    ullage_frac::Float64
    innertank() = new()
end

"""
doublewalled_tank calcualtes the weight and heat transfer of a 
double walled LH2 tank pg 375 of Baron
## Inputs
fuse_clearance: [m] Clearance from fuselage
d_vaccuum     : [m] distance between inner and outer vessel
Ltank         : [m] Outer most length of tank
Δp             : [Pa] Internal storage pressure
σa_inner      : [Pa] Allowable stress in inner vessel ~129.2 MPa for SS 316
ew            : [-]  Weld efficiency - Fully radiographed -> 1.0 (See Barron pg 361)
ARtank            : [-]  Aspect ratio of heads - 1 => hemisphere
ρinner        : [kg/m³] Density of inner tank material. High Nickel steel ideal to avoid embrittilement (SS314 = 7900 kg/m³)
Nstiff_in     : [-] Number of inner vessel stiffners
θin_support   : [rad] Location of inner vessel support 
σa_outer      : [Pa] Allowable stress in outer vessel
poiss         : [-] Poission ratio of outer vessel material
E_outer       : [Pa] Young's modulus of outer vessel material
ρouter        : [kg/m³] Density of outer tank material. 
ullage        : [-]   Ullage fraction - empty space in the tank for vapor
ρfuel         : [kg/m³] fuel density
subscript 1 represent inner tank, 2 is outer tank
use of inner or outer for R (radius) is the inner and outer Radius of a wall
"""
function doublewalled_tank( Rfuse::Float64, dRfuse::Float64, clearance_fuse::Float64, 
            Wfuel::Float64, ρfuel::Float64, ARtank::Float64, ptank::Float64, ullage_frac::Float64,
            siginner::Float64, ρinner::Float64,ew::Float64,
            material_insul::Vector{String}, t_cond::Array{Float64,1}, 
            m_boiloff::Float64,
            σa_outer::Float64, ρouter::Float64, E_outer::Float64,
            Nstiff_in, θin_support, θout_support1,θout_support2, poiss, ftankadd, ullage_frac)
  
    
    #Create structure with inner tank parameters for tankWmech()
    inner_tank = innertank() #instantiate tank object
    inner_tank.ftankstiff = 0.0 #stiffeners are sized explicitly later
    inner_tank.ftankadd = ftankadd
    inner_tank.ptank = ptank
    inner_tank.sigskin = siginner
    inner_tank.material_insul = material_insul
    inner_tank.rhoskintank = ρinner
    inner_tank.clearance_fuse = clearance_fuse
    inner_tank.ARtank = ARtank
    inner_tank.ullage_frac = ullage_frac

    _, l_cyl1, t_cyl1, Rtank1_outer, Vfuel, Wtank1, _, Winsul_sum, 
    t_head1, Whead1, Wcyl1, Winsul, S_inner, Shead_insul, _ =
                tankWmech(inner_tank, t_cond, ρfuel,
                Rfuse, dRfuse, wfb, nfweb,
                Wfuel)

    
## Stiffening
    ϕlist, klist, ϕmax, kmax = stiffeners_bendingM(θin_support) #TODO check this function
    Mmax = kmax * Rtank1_outer* (W_inner/Nstiff_in) /(2π) #use outer dia for now
    Zmax = Mmax/σa_inner
        
    #Use standard metric beams to get weights per unit length #TODO these discontinuities may cause numerical issues
    if Zmax < 89.9e-6
        # W 100 x 100 x 19.3
        W′ = 19.3 # kg/m
    elseif Zmax < 139.5e-6
        # W 130 x 130 x 23.8
        W′ = 23.8 # kg/m
    else
        W′=27.5 # 13 x 13 x 27.5
    end
       
    W_stiffners1 = Nstiff_in * 2π*Rtank1_outer * W′ * gee
    Wtank1 = W_stiffners1 + (Wcyl1 + Whead1) # update with stiffener weight

    W_inner = Wfuel_tot + Wtank1

    ############### Outer vessel ############### 

    # Collapsing pressure for "long" cylinder
    pc = 4*pSL
    
    Rtank2_outer = Rfuse - fuse_clearance 
    N_outer_rings= optimize_outer(l_cyl1,Rtank2_outer,poiss,ρouter, σa_outer,E_outer,pc,W_inner,θout_support1,θout_support2,ARtank)

    Wtank2,W_stiffners2,Wcyl2,
    Whead2,S_outer,Shead2,Scyl2,
    t_cyl2,t_head2 = size_outer_tank(N_outer_rings,l_cyl1,Rtank2_outer,poiss,ρouter,σa_outer, E_outer,pc,W_inner,θout_support1,θout_support2,ARtank)

    Rtank2_inner = Rtank2_outer-t_cyl2

    l_tank = l_cyl1+ 2*( (Rfuse-fuse_clearance) /ARtank)

    Wtank = Wtank1 + Wtank2 + Winsul_sum 
    # println(" Wtank1 $Wtank1 Wtank2 $Wtank2 W_stiffners1 $W_stiffners1 W_stiffners2 $W_stiffners2 ")

    Wtank_total = Wtank + Wfuel_tot
    η = Wfuel_tot/ Wtank_total *100

    # println("Wtank2 = $Wtank2")
    #     println("Tank total = $Wtank_total, η = $η")
    # println("fadd = ", W_stiffners1/(Wtank_total - W_stiffners1))
    # println((Wtank_total - W_stiffners1)*(1+0.032413))

    return Wtank_total, Wtank, Vfuel, η, 
            Whead1, Wcyl1, Wtank1, W_stiffners1, 
            Whead2, Wcyl2, Wtank2, W_stiffners2, N_outer_rings, Winsul_sum,
            l_cyl1, l_tank, t_cyl1, t_head1, t_cyl2,t_head2, 
            Rtank1_outer, Rtank1_inner,  Rtank2_outer, Rtank2_inner,
            S_inner,S_outer, Shead_insul, Shead1, Scyl1, Shead2, Scyl2

end

""" find min of outer tank
and return  Number of stiffeners and relevant outer tank params from size_outer_tank
"""
function optimize_outer(l_cyl1,Rtank2_outer,poiss,ρouter,σa_outer, E_outer,pc,W_inner,θout_support1,θout_support2,ARtank)
    N0=round(l_cyl1/2)
    W0,_ = size_outer_tank(N0,l_cyl1,Rtank2_outer,poiss,ρouter, σa_outer,E_outer,pc,W_inner,θout_support1,θout_support2,ARtank)
    # println("startin inner", N0)
    Nopt=N0 # guess
    for N =2:round(l_cyl1/0.5)
        Wnew,_=size_outer_tank(N,l_cyl1,Rtank2_outer,poiss,ρouter,σa_outer, E_outer,pc,W_inner,θout_support1,θout_support2,ARtank)
        # println("W0, $W0, Wnew, $Wnew  $N")
        if Wnew<W0
            # minimum found
            # println("yeeeeeeee  ",N)
            Nopt=N
            W0=Wnew
        end
    end
    # if cant find a result
    return Nopt
end

function size_outer_tank(N,l_cyl1,Rtank2_outer,poiss,ρouter,σa_outer, E_outer,pc,W_inner,θout_support1,θout_support2,ARtank)
    #N number of stiffeners
    L = l_cyl1/(N-1)
    Doo=2*Rtank2_outer
    t_cyl2 = Doo* ( pc*(1-poiss^2)^0.75 *L/Doo/(2.42*E_outer)  )^(2/5)    # guess using long cylinder

    if L/Doo < 1.140*(1-poiss^2)^(1.0/4)*sqrt(Doo/t_cyl2)  # Short cylinder (see barrons pg366 eq7.8)
        # println("Warning: Outer Cycl isn't really 'long' - using 'short' tank relation instead")
        t_Dguess=t_cyl2/Doo
        t_cyl2 = tank_buckling(pc, E_outer, poiss, L, Doo, t_Dguess)
        # println("t_cyl2 = $t_cyl2")
    else # Long Cylinder
        t_cyl2 = Doo*  (  pc*(1-poiss^2)^0.75 *L/Doo/(2.42*E_outer)     )^(2/5)    #
    end

    if ARtank==2.0
        K1=0.90# See table 7.6 for D/D1=2.0 in Barron p. 367
    elseif ARtank==1.0
        K1=0.50
    else  
        println("ARtank of heads not supported, see doublewalltank")
        K1=1.0
    end
    t_head2 = K1*Doo*sqrt(pc*sqrt(3*(1 - poiss^2))/ (0.5*E_outer)  )

    ## Areas
    Shead2 = 2π*Rtank2_outer^2 * (1.0/3.0 + 2.0/3.0*(1/ARtank)^1.6) ^ (1/1.6) # Surface area of head
    Scyl2  = 2π*Rtank2_outer*l_cyl1  # Cross-sectional area

    S_outer=Shead2+Scyl2 

    ## Volume and Weight
    Vcyl  = Scyl2*t_cyl2
    Vhead = Shead2*t_head2

    Wcyl2  = Vcyl *ρouter*gee
    Whead2 =  2*Vhead*ρouter*gee


    ϕlist, klist, ϕmax, kmax_outer = stiffeners_bendingM_outer(θout_support1,θout_support2)
    Mmax_outer = kmax_outer * Rtank2_outer* (W_inner/N) /(2π) #use outer dia for now
    Zmax = Mmax_outer/σa_outer
    
    #Use standard metric beams to get weights per unit length
    if Zmax < 89.9e-6
        # W 100 x 100 x 19.3
        W′ = 19.3 # kg/m
    elseif Zmax < 139.5e-6
        # W 130 x 130 x 23.8
        W′ = 23.8 # kg/m
    else
        W′=27.5 # 13 x 13 x 27.5
    end
    # println("W' $(W′)")
    
    W_stiffners2 = N * 2π*Rtank2_outer * W′ * gee

    Wtank2 = Wcyl2 + Whead2 + W_stiffners2

    return Wtank2,W_stiffners2,Wcyl2,Whead2,S_outer,Shead2,Scyl2,t_cyl2,t_head2
end



function tank_buckling(pc, E, poiss, L, Do, t_Dguess)
    f(t_D) = 2.42*E*(t_D)^(5.0/2) / ( (1-poiss^2)^(3.0/4.0) * ((L/Do) - 0.45*sqrt(t_D)) ) - pc
    try
        t_D = find_zero(f, t_Dguess)
        return t = t_D*Do
    catch
        println("fail t find t_cyl2 using short cylinder method")
        return t_Dguess
    end
    
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
            k[i] = 0.5*cos(ϕ) - (π - ϕ)*sin(ϕ) + θ  + cos(θ) + cos(ϕ)*(sin(θ)^2)
        end
    end
    kmax, imax = findmax(abs.(k))
    ϕmax = ϕlist[imax]
    return ϕlist, k, ϕmax, kmax
end

function stiffeners_bendingM_outer(θ1,θ2)
    ϕlist = LinRange(0.0, π, 1800)
    k = zeros(length(ϕlist))

    for (i,ϕ) in enumerate(ϕlist)
        if 0 ≤ ϕ ≤ θ1
            k[i] = cos(ϕ)*(sin(θ2)^2 - sin(θ1)^2 ) + (cos(θ2) - cos(θ1) ) 
                         - (π-θ2)*sin(θ2) +  (π-θ1)*sin(θ1)
        elseif θ1 ≤ ϕ ≤ θ2
            k[i] = cos(ϕ)*(sin(θ2)^2 - sin(θ1)^2 ) + (cos(θ2) - cos(θ1) ) 
                         - (π-θ2)*sin(θ2) +  π*sin(ϕ) -θ1*sin(θ1)
        elseif θ2 ≤ ϕ ≤ π
            k[i] = cos(ϕ)*(sin(θ2)^2 - sin(θ1)^2 ) + (cos(θ2) - cos(θ1) )  +
                          (θ2*sin(θ2) - θ1*sin(θ1)  )
        end
    end
    kmax, imax = findmax(abs.(k))
    ϕmax = ϕlist[imax]
    return ϕlist, k, ϕmax, kmax
end



function tank_thermal_vacumm(l_cyl::Float64  , Rtank1_outer::Float64,Rtank2_inner::Float64, 
    hconvgas::Float64, h_LH2::Float64,  hconvair::Float64, h_v::Float64,
    t_cond::Array{Float64,1}, k::Array{Float64,1},
    Tfuel::Float64 , Tair::Float64, time_flight::Float64, 
    S_outer::Float64, S_inner::Float64, Sheads_insul::Array{Float64,1}, p_vaccum::Float64,
     a_outer::Float64, a_inner::Float64 )

    N = length(t_cond)       # Number of layers in insulation
    thickness = sum(t_cond)  # total thickness of insulation

    #--- Heat flux and resistances
    ΔT = Tair - Tfuel  # Overall temperature drop between ambient and LH2

    qfac = 1.3         # Account for heat leak from pipes and valves

    # [TODO] Move constants to constant file
    σ = 5.67e-8
    ε = 0.95    # white aircraft (Verstraete)

    hradair = σ * ε * ((Tair^2) + (Tfuel^2)) * (Tair + Tfuel) #radiative coeff
    h_air = hconvair + hradair # Combines radiative and convective heat transfer at outer end
    Rair_conv_rad = 1 / (h_air * S_outer  )  # thermal resistance of ambient air (incl. conv and rad)
    r_inner = Rtank1_outer #- thickness

    R_LH2 = 1 / (h_LH2 * S_inner) #thermal resistance of LH2

    R_mli      = zeros(Float64, N)  #size of foam and vapor barrier 
    R_mli_ends = zeros(Float64, N)
    R_mli_cyl  = zeros(Float64, N)

    for i in 1:N
        if k[i]==0 # if vacuum tank skip it
            break
        end
        R_mli_cyl[i]  = log((r_inner  + t_cond[i])/ (r_inner)) / (2π*l_cyl * k[i]) #Resistance of each MLI layer
        R_mli_ends[i] = t_cond[i] / (k[i] * 2*Sheads_insul[i])
        # Parallel addition of resistance
        R_mli[i]  = (R_mli_ends[i] * R_mli_cyl[i]/(R_mli_ends[i] + R_mli_cyl[i])) 

        # Update r_inner
        r_inner   = r_inner + t_cond[i]  
    end

    R_mli_tot = sum(R_mli)  #Total thermal resistance of foam and vapor barrier

    #### If vac tank need to account for imperfect vacuum 
    mu=18.47*10^-6 #mu at 300K
    gc= 1.0 # molecular weight of air g/mole 
    Rgas=287.05  # specific gas constant
    lambda=mu/p_vaccum*sqrt(pi*Rgas*Tair/2/gc)
    gamma=1.4

    vac_index=N-1
    if lambda>Rtank2_inner-Rtank1_outer
        Fa=1/( 1/a_inner+(S_outer/S_inner)*(1/a_outer-1)  )
        G=(gamma+1)/(gamma-1)*sqrt(gc*Rgas/(8*pi*Tfuel))*Fa
        R_vac_conv = 1/(G*p_vaccum * S_inner  )  # convective resistance due to imperfect vac
        R_vac_rad = 1/(hradair * S_inner  )  # radiative resistance 
        R_vac_sum = R_mli[vac_index]*R_vac_conv*R_vac_rad/(R_mli[vac_index]+R_vac_conv+R_vac_rad) # parallel sum of MLI, conv and raditation
        # println("$R_vac_sum  $Rair_conv_rad $( R_vac_conv * R_vac_rad/( R_vac_conv + R_vac_rad )) ")

    else
        println("no conv")
        R_vac_conv=0.0
        R_vac_rad = 1 / (hradair * S_inner  ) # radiative resistance
        R_vac_sum = R_mli[vac_index]*R_vac_rad/(R_mli[vac_index]+R_vac_rad) # parallel sum of MLI, conv and raditation
    end

    Req = Rair_conv_rad + R_vac_sum+ R_mli_tot + R_LH2   # Total equivalent resistance of thermal circuit
    # println("$Rair_conv_rad +$R_vac_sum + $R_mli_tot  + $R_LH2")

    q = qfac * ΔT / Req     # Heat flux from ambient to LH2, including extra heat leak from valves etc as in eq 3.20 by Verstraete
    mdot_boiloff = q / h_v  # Boil-off rate equals the heat flux divided by heat of vaporization
    m_boiloff = mdot_boiloff * time_flight # Boil-off mass calculation

    return  m_boiloff, mdot_boiloff, Req
end





function tanksizeVac(Rfuse::Float64, dRfuse::Float64, fuse_clearance::Float64, 
                    Wfuel::Float64, ρfuel::Float64, ARtank::Float64, Δp::Float64, ullage::Float64,
                    time_flight::Float64, threshold_percent::Float64,
                    σa_inner::Float64, ρinner::Float64,ew::Float64,
                    Tair::Float64, Tfuel::Float64, 
                    h_LH2::Float64, h_v::Float64, hconvair::Float64, hconvgas::Float64,
                    p_vaccum_i::Float64, 
                    σa_outer::Float64, ρouter::Float64, E_outer::Float64,
                    rho_insul::Array{Float64,1}, t_cond::Array{Float64,1},  k::Array{Float64,1}, 
                    a_outer::Float64,a_inner::Float64, Nstiff_in::Float64, 
                    θin_support::Float64, θout_support1::Float64,θout_support2::Float64, poiss::Float64)

    Wfuel_init = Wfuel
    # not accurate
    m_boiloff = 0. #initial value of boil-off mass
    Req=0.
    mdot_boiloff = 0.0
    p_vaccum=p_vaccum_i

    Wtank_total,Wfuel_tot, Wtank, Vfuel, η, 
    Whead1, Wcyl1, Wtank1, W_stiffners1, 
    Whead2, Wcyl2, Wtank2, W_stiffners2, N_outer_rings, Winsul_sum,
    l_cyl1, L_tank, t_cyl1, t_head1, t_cyl2,t_head2, 
    Rtank1_outer, Rtank1_inner,  Rtank2_outer, Rtank2_inner,
    S_inner,S_outer, Sheads_insul, Shead2,Scyl2,Shead1,Scyl1 = doublewalled_tank( Rfuse, dRfuse, fuse_clearance, 
                                                                Wfuel, ρfuel, ARtank, Δp, ullage,
                                                                σa_inner, ρinner,ew,
                                                                rho_insul, t_cond, 
                                                                m_boiloff,
                                                                σa_outer, ρouter, E_outer,
                                                                Nstiff_in, θin_support, θout_support1,θout_support2, poiss)
    #optimize boil off mass according to threshold
    while p_vaccum >10^-5  # up to Very High Order Vac
        m_boiloff, mdot_boiloff,Req = tank_thermal_vacumm(l_cyl1  ,Rtank1_outer,Rtank2_inner, 
                                                    hconvgas, h_LH2,  hconvair, h_v,
                                                    t_cond, k,
                                                    Tfuel , Tair, time_flight, 
                                                    S_outer, S_inner,Sheads_insul, p_vaccum,
                                                    a_outer , a_inner )


        Wtank_total,Wfuel_tot, Wtank, Vfuel, η, 
        Whead1, Wcyl1, Wtank1, W_stiffners1, 
        Whead2, Wcyl2, Wtank2, W_stiffners2, N_outer_rings, Winsul_sum,
        l_cyl1, L_tank, t_cyl1, t_head1, t_cyl2,t_head2, 
        Rtank1_outer, Rtank1_inner,  Rtank2_outer, Rtank2_inner,
        S_inner,S_outer, Sheads_insul, Shead2,Scyl2,Shead1,Scyl1 =doublewalled_tank( Rfuse, dRfuse, fuse_clearance, 
                                                    Wfuel, ρfuel, ARtank, Δp, ullage,
                                                    σa_inner, ρinner,ew,
                                                    rho_insul, t_cond, 
                                                    m_boiloff,
                                                    σa_outer, ρouter, E_outer,
                                                    Nstiff_in, θin_support, θout_support1,θout_support2, poiss)
        Wfuel = Wfuel_init
        if(mdot_boiloff*gee*3600/Wfuel *100) < threshold_percent
            # println("iter = $n")
            # println("Boil off fraction $(mdot_boiloff*gee*3600/Wfuel *100) %")
            break
        end
        p_vaccum=p_vaccum*0.5
        # println(p_vaccum)
        
    end

    if(mdot_boiloff*gee*3600/Wfuel *100) > threshold_percent
        # println("iter = $n")
        println("$p_vaccum Pa cant meet boiloff Vac tank")
    end
    return Wtank_total,η, l_cyl1, S_inner,S_outer, 
    Shead2,Scyl2,Shead1,Scyl1,Wtank1,
    Vfuel, Wtank, N_outer_rings, W_stiffners2,Wcyl2,
    Whead2, Winsul_sum,Wfuel_tot, t_cyl1,t_head1, 
    t_cyl2,t_head2,  L_tank, Rtank2_outer, Rtank2_inner,Rtank1_outer, Rtank1_inner,
    m_boiloff, mdot_boiloff ,p_vaccum,Req

end

const pSL = 101325

ρfuel = 70.0
Rfuse = 2.5
dRfuse = 0.3
Tfuel = 20.0
time_flight = 7*3600.0
wfb = 0.0
nfweb = 1.0
Wfuel = 1e5
ifuel = 40
z = 11e3
Mair = 0.8
xftank = 15.0
hconvgas = 0.0
Tair  = 288.0 #Heated cabin temp
h_v = 447000.0
h_LH2 = 210.0
hconvair = 15.0
threshold_percent  = 0.1

σa_inner = 172.4e6 ;  ρinner =  2825.0  #AL 2219 Brewer / energies stress for operating conditions (290e6 ultimate operatoin)

# println("Using Vacumm tank")
poiss = 0.33

k_MLI =  0.00015# 
rho_MLI = 3.35 #kg/m2  Table 9 A139 40 layers 0.94 layer/mm http://etd.fcla.edu/CF/CFE0003419/Johnson_Wesley_L_201012_MS.pdf 
# rho_MLI = 

σa_outer = 172.4e6  ; E_outer = 73.8e9; ρouter = 2825.0
# σa_outer = 332.0e6; E_outer = 164.0e9; ρouter = 7870.0

Nstiff_in = 3.0
θin_support = 80*π/180.0
θout_support1 = 70*π/180.0 ; θout_support2 = 110*π/180.0
a_outer=0.29; a_inner=0.59 # accomodation coefficient of Helium/Hydrogen 

p_vaccum_i = 0.01 #Pa #initial vac level
d_vaccuum=0.065 # vacumm size recommended by brewer (brewer book 191/204)

# closed cell foam, vapor bar, Vac/MLI, vapor barr
t_cond = [ 0.015, 1.524e-5,     d_vaccuum, 1.524e-5] #brewer modified
k = [5e-3,      5e-3,            k_MLI,  5e-3,  ] #foam conductivities
rho_insul = [35.24, 14764,  rho_MLI,14764] #energies

ew = 0.9 #weld efficiency
ARtank = 2.0
fuse_clearance = 0.10
ρfuel = parg[igrhofuel] 
ptank = 2.0*101325.0 #atm
ftankstiff = 0.1
ftankadd   = 0.1
ullage = 0.1


Wtank_total,η, l_cyl1, S_inner,S_outer, 
Shead2,Scyl2,Shead1,Scyl1,Wtank1,
Vfuel, Wtank, N_outer_rings, W_stiffners2,Wcyl2,
Whead2, Winsul_sum ,Wfuel_tot, t_cyl1,t_head1, 
t_cyl2,t_head2,  L_tank, Rtank2_outer, Rtank2_inner,Rtank1_outer, Rtank1_inner,
m_boiloff, mdot_boiloff ,p_vaccum,Req = tanksizeVac(Rfuse, dRfuse, fuse_clearance, 
                                            Wfuel, ρfuel, ARtank, ptank, ullage,
                                            time_flight, threshold_percent,
                                            σa_inner, ρinner,ew,
                                            Tair, Tfuel, 
                                            h_LH2, h_v, hconvair, hconvgas,
                                            p_vaccum_i,
                                            σa_outer, ρouter, E_outer,
                                            rho_insul, t_cond, k,
                                            a_outer,a_inner,Nstiff_in, 
                                            θin_support, θout_support1,θout_support2, poiss)