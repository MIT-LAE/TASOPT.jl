using Roots
using NLopt

"""
    stiffeners_bendingM(θ)
This function can be used to calculate the bending moment distribution in a stiffener ring for an inner cryogenic tank.
It applies Eqs. (7.4) and (7.5) in Barron (1985) to find the bending moment distribution. The function returns the
maximum value of ``k = 2πM/(WR)`` on the ring's circumference.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `θ::Float64`: angular position of tank supports, measured from the bottom of the tank (rad).

    **Outputs:**
    - `ϕmax::Float64`: angular position of maximum bending moment on ring circumference (rad).
    - `kmax::Float64`: Maximum value of the ratio ``k = 2πM/(WR)`` on ring circumference.
"""
function stiffeners_bendingM(θ)
    ϕlist = LinRange(0.0, π, 180)
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
    return ϕmax, kmax
end

"""
    stiffeners_bendingM_outer(θ1,θ2)
This function can be used to calculate the bending moment distribution in a stiffener ring for an outer cryogenic tank.
It applies Eqs. (7.13)-(7.15) in Barron (1985) to find the bending moment distribution. The function returns the
maximum value of ``k = 2πM/(WR)`` on the ring's circumference.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `θ1::Float64`: angular position of bottom tank supports, measured from the bottom of the tank (rad).
    - `θ2::Float64`: angular position of top tank supports, measured from the bottom of the tank (rad).

    **Outputs:**
    - `ϕmax::Float64`: angular position of maximum bending moment on ring circumference (rad).
    - `kmax::Float64`: Maximum value of the ratio ``k = 2πM/(WR)`` on ring circumference.
"""
function stiffeners_bendingM_outer(θ1,θ2)
    ϕlist = LinRange(0.0, π, 180)
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
    return ϕmax, kmax
end

"""
    stiffener_weight(tanktype, W, Rtank, s_a, ρstiff, θ1, θ2 = 0.0, Nstiff = 2.0, l_cyl = 0, E = 0)
This function calculates the weight of a single stiffener in an inner or outer tank for a given inner tank weight.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `tanktype::String`: type of tank, options are "inner" or "outer".
    - `W::Float64`: load carried by a stiffener ring (N).
    - `Rtank::Float64`: tank radius (m).
    - `s_a::Float64`: maximum allowable stress in stiffener material (Pa).
    - `ρstiff::Float64`: stiffener density (kg/m^3).
    - `θ1::Float64`: angular position of bottom tank supports, measured from the bottom of the tank (rad).
    - `θ2::Float64`: angular position of top tank supports, measured from the bottom of the tank (rad). Only used with "outer" tank.
    - `Nstiff::Float64`: total number of stiffeners on outer tank. Only used with "outer" tank.
    - `l_cyl::Float64`: length of cylindrical portion of tank (m). Only used with "outer" tank.
    - `E::Float64`: Young's modulus of stiffener material (Pa). Only used with "outer" tank.

    **Outputs:**
    - `Wstiff::Float64`: weight of a stiffener ring (N).
"""
function stiffener_weight(tanktype, W, Rtank, s_a, ρstiff, θ1, θ2 = 0.0, Nstiff = 2.0, l_cyl = 0, E = 0)
    
    if tanktype == "inner" 
        _, kmax = stiffeners_bendingM(θ1) #Find k = 2πM/(WR)
        Icollapse = 0 #Inner tank cannot collapse as it is pressurized

    elseif tanktype == "outer"
        _, kmax = stiffeners_bendingM_outer(θ1, θ2) #Find k = 2πM/(WR)
        pc = 4 * pref #Critical pressure is 4 times atmospheric pressure, Eq. (7.11) in Barron (1985)
        Do = 2 * Rtank #outer diameter

        L = l_cyl/ (Nstiff - 1) #Length of portion between supports
        Icollapse = pc * Do * L / (24 * E) #Second moment of area needed to avoid collapse
    end

    Mmax = kmax * W * Rtank / (2π) #Maximum bending moment due to loads
    Z = Mmax / s_a #required section modulus to withstand bending loads

    #Assume sectional properties of a 100 x 100 I-beam
    W = 100e-3 #Flange width
    t_w = 7.1e-3 #Web thickness
    t_f = 8.8e-3 #Flange thickness

    #For an I-beam, I > W * H^2 * t_f / 2
    #The required second moment of area is I = Icollapse + Z * H/2
    #Find beam height by solving W * H^2 * t_f / 2 = Icollapse + Z * H/2

    a = t_f * W / 2 #Coefficients in quadratic equation
    b = -Z/2
    c = -1 * Icollapse

    H = (-b + sqrt(b^2 - 4 * a * c)) / (2 * a) #Solve quadratic eq.
    S = 2 * W * t_f + (H - t_w) * t_f #Beam cross-sectional area

    Wstiff = gee * ρstiff * S * 2π * Rtank #Weight of a single stiffener running along circumference
    return Wstiff
end

"""
    optimize_outer_tank(fuse_tank, Winnertank, l_cyl, θ1, θ2)
This function optimizes the number of intermediate stiffener rings to minimize the weight of an outer tank.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuse_tank::Struct`: structure with tank parameters.
    - `Winnertank::Float64`: weight of inner tank and contents (N).
    - `l_cyl::Float64`: length of cylindrical portion of tank (m).

    **Outputs:**
    - `Ninterm::Float64`: optimum number of intermediate stiffener rings.
"""
function optimize_outer_tank(fuse_tank, Winnertank, l_cyl)

    obj(x, grad) = size_outer_tank(fuse_tank, Winnertank, l_cyl, x[1])[1] #Minimize Wtank
    
    initial_x = [1.0]
    
    #Set bounds
    lower = [0.0]
    upper = [50.0]
   
    #Use NLopt.jl to minimize function 
    opt = Opt(:LN_NELDERMEAD, length(initial_x))
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.ftol_rel = 1e-9
    opt.maxeval = 100  # Set the maximum number of function evaluations

    opt.min_objective = obj

    (minf,xopt,ret) = NLopt.optimize(opt, initial_x)
    Ninterm = xopt[1]
    return Ninterm
end


mutable struct outertank
    poissouter::Float64
    Eouter::Float64
    UTSouter::Float64
    rhoouter::Float64
    ftankadd::Float64
    wfb::Float64
    nfweb::Float64
    ARtank::Float64
    Rfuse::Float64
    clearance_fuse::Float64
    theta_outer::Vector{Float64}
    outertank() = new()
end
 
const pref = 101325
const gee = 9.81

l_cyl = 5

fuse_tank = outertank()
fuse_tank.poissouter = 0.3
fuse_tank.Eouter = 73e9
fuse_tank.rhoouter = 2840
fuse_tank.UTSouter = 470e6
fuse_tank.ftankadd = 0.1
fuse_tank.wfb = 0.0
fuse_tank.nfweb = 1.0
fuse_tank.ARtank = 2.0
fuse_tank.Rfuse = 2.2
fuse_tank.clearance_fuse = 0.1
fuse_tank.theta_outer = [1.0, 2.0]

"""
    size_outer_tank(fuse_tank, Winnertank, l_cyl, Ninterm)
This function sizes the outer tank and calculates the weights of its components.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuse_tank::Struct`: structure with tank parameters.
    - `Winnertank::Float64`: weight of inner tank and contents (N).
    - `l_cyl::Float64`: length of cylindrical portion of tank (m).
    - `Ninterm::Float64`: optimum number of intermediate stiffener rings.

    **Outputs:**
    - `Wtank::Float64`: total weight of outer tank (N).
    - `Wcyl::Float64`: weight of cylindrical portion of outer tank (N).
    - `Whead::Float64`: weight of one elliptical outer-tank head (N).
    - `Wstiff::Float64`: total weight of stiffener material (N).
    - `S_outer::Float64`: surface area of outer tank (m^2).
    - `Shead::Float64`: surface area of one outer tank head (m^2).
    - `Scyl::Float64`: surface area of cylindrical portion of tank (m^2).
    - `t_cyl::Float64`: wall thickness of cylindrical portion of tank (m).
    - `t_head::Float64`: wall thickness of tank head (m). 
"""
function size_outer_tank(fuse_tank, Winnertank, l_cyl, Ninterm)
    #Unpack parameters in fuse_tank
    poiss = fuse_tank.poissouter
    Eouter = fuse_tank.Eouter
    ρouter = fuse_tank.rhoouter
    UTSouter = fuse_tank.UTSouter
    ftankadd = fuse_tank.ftankadd
    wfb = fuse_tank.wfb
    nfweb = fuse_tank.nfweb
    ARtank = fuse_tank.ARtank
    θ_outer = fuse_tank.theta_outer

    θ1 = θ_outer[1]
    θ2 = θ_outer[2]

    Nmain = 2 #Tanks typically have two main support rings
    pc = 4 * pref #4*p_atm; Collapsing pressure, Eq. (7.11) in Barron (1985)
    s_a = UTSouter / 4

    #Calculate outer tank geometry
    Rtank_outer = fuse_tank.Rfuse - fuse_tank.clearance_fuse
    Do = 2 * Rtank_outer #outside diameter

    Nstiff = Nmain + Ninterm #Total number of stiffeners
    L = l_cyl / (Nstiff - 1) #There are two stiffeners at the ends, so effective number of sections in skin is N - 1
    L_Do = L / Do

    #Find cylinder wall thickness. This applies to a short cylinder.
    pressure_res(t_D) = 2.42*Eouter*(t_D)^(5/2) / ( (1 - poiss^2)^(3/4) * (L_Do - 0.45*sqrt(t_D)) ) - pc
    t_Do = find_zero(pressure_res, 1e-3) #Find root with Roots.jl
    t_cyl = t_Do * Do

    #Find head wall thickness
    if ARtank == 2.0
          K1 = 0.90# See table 7.6 for D/D1=2.0 in Barron p. 367
    elseif ARtank == 1.0
          K1 = 0.50
    else  
          println("ARtank of heads not supported, see size_outer_tank()")
          K1=1.0
    end
    t_head = K1 * Do * sqrt(pc * sqrt(3*(1 - poiss^2))/ (0.5*Eouter))

    ## Areas
    wfblim = max( min( wfb , Rtank_outer) , 0.0 )
    thetafb = asin(wfblim / Rtank_outer)

    Shead = (2.0*π + 4.0*nfweb*thetafb)*(Rtank_outer)^2* (0.333 + 0.667*(1/ARtank)^1.6 )^0.625
    Scyl  = 2π*Rtank_outer*l_cyl  # Cross-sectional area

    Souter = Scyl + 2*Shead

    ## Volume and Weight
    Vcyl  = Scyl*t_cyl
    Vhead = Shead*t_head

    Wcyl  = Vcyl*ρouter*gee
    Whead =  Vhead*ρouter*gee

    Wtank_no_stiff = Wcyl + 2 * Whead

    # Size stiffeners
    tanktype = "outer"

    Wmainstiff = stiffener_weight(tanktype, Winnertank / Nmain, Rtank_outer, #Weight of one main stiffener, each one 
                                s_a, ρouter, θ1, θ2, Nstiff, l_cyl, Eouter)  #carries half the inner tank load
                                                                            
    Wintermstiff = stiffener_weight(tanktype, 0.0, Rtank_outer, 
                                s_a, ρouter, θ1, θ2, Nstiff, l_cyl, Eouter) #Weight of one intermediate stiffener, which carries no load

    Wstiff = Nmain * Wmainstiff + Ninterm * Wintermstiff #Total stiffener weight

    Wtank = (Wtank_no_stiff + Wstiff) * (1 + ftankadd) #Find total tank weight, including additional mass factor

    return Wtank, Wcyl, Whead, Wstiff, Souter, Shead, Scyl, t_cyl, t_head
end


