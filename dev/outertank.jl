using Roots
using NLopt

"""
## Stiffening rings
returns 2πM/(WR)
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

function optimize_inner_stiffeners(tanktype, Wtanktotal, Rtank, s_a, ρstiff)
    obj(x, grad) = stiffener_weight(x, tanktype, Wtanktotal, Rtank, s_a, ρstiff, l_cyl, E) #Opjective function is total stiffener weight
    
    if tanktype == "inner"
        initial_x = [2.0, 1.0]
        #Set bounds
        lower = [2.0, 0.0]
        upper = [50.0, pi/2]
    elseif tanktype == "outer"
        initial_x = [2.0, 1.0, 2.0]
        #Set bounds
        lower = [2.0, 0.0, pi/2]
        upper = [50.0, pi/2, pi]
    end
    
    #Use NLopt.jl to minimize function 
    opt = Opt(:LN_NELDERMEAD, length(initial_x))
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.ftol_rel = 1e-5
    opt.maxeval = 100  # Set the maximum number of function evaluations

    opt.min_objective = obj
    println(obj(initial_x, 0))

    (minf,xopt,ret) = NLopt.optimize(opt, initial_x)
    return xopt
end

function stiffener_weight(tanktype, W, Rtank, s_a, ρstiff, Nstiff, θ1, θ2 = 0.0, l_cyl = 0, E = 0)
    #Unpack optimization variables
    if tanktype == "inner"
        _, kmax = stiffeners_bendingM(θ1)
        Icollapse = 0
    elseif tanktype == "outer"
        _, kmax = stiffeners_bendingM_outer(θ1, θ2)
        pc = 4 * pref
        Do = 2 * Rtank

        L = l_cyl/ (Nstiff - 1)
        Icollapse = pc * Do * L / (24 * E)

    end

    Mmax = kmax * W * Rtank / (2π)
    Z = Mmax / s_a #required section modulus

    W = 100e-3
    t_w = 7.1e-3
    t_f = 8.8e-3

    #For an I-beam, I > W * H^2 * t_f / 2
    a = t_f * W / 2 
    b = -Z/2
    c = -1 * Icollapse

    H = (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
    S = 2 * W * t_f + (H - t_w) * t_f

    Wstiff = gee * ρstiff * S * 2π * Rtank
    return Wstiff
end

function optimize_outer_tank(fuse_tank, Winnertank, l_cyl, θ1, θ2)
    obj(x, grad) = size_outer_tank(fuse_tank, Winnertank, l_cyl, x[1], θ1, θ2)[1] #Minimize Wtank
    
    initial_x = [1.0]
    #Set bounds
    lower = [0.0]
    upper = [50.0]
   
    #Use NLopt.jl to minimize function 
    opt = Opt(:LN_NELDERMEAD, length(initial_x))
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.ftol_rel = 1e-5
    opt.maxeval = 100  # Set the maximum number of function evaluations

    opt.min_objective = obj

    (minf,xopt,ret) = NLopt.optimize(opt, initial_x)
    return xopt
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
fuse_tank.Rfuse = 2.5
fuse_tank.clearance_fuse = 0.1

function size_outer_tank(fuse_tank, Winnertank, l_cyl, Ninterm, θ1, θ2)
    #Unpack parameters in fuse_tank
    poiss = fuse_tank.poissouter
    Eouter = fuse_tank.Eouter
    ρouter = fuse_tank.rhoouter
    UTSouter = fuse_tank.UTSouter
    ftankadd = fuse_tank.ftankadd
    wfb = fuse_tank.wfb
    nfweb = fuse_tank.nfweb
    ARtank = fuse_tank.ARtank

    Nmain = 2 #Tanks typically have two main support rings
    pc = 4 * pref #4*p_atm; Collapsing pressure, Eq. (7.11) in Barron (1985)
    s_a = UTSouter / 4

    #Calculate outer tank geometry
    Rtank_outer = fuse_tank.Rfuse - fuse_tank.clearance_fuse
    Do = 2 * Rtank_outer #outside diameter

    L = l_cyl / (Nmain + Ninterm - 1) #There are two stiffeners at the ends, so effective number of sections in skin is N - 1
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

    S_outer = Shead + 2 * Scyl

    ## Volume and Weight
    Vcyl  = Scyl*t_cyl
    Vhead = Shead*t_head

    Wcyl  = Vcyl*ρouter*gee
    Whead =  Vhead*ρouter*gee

    Wtank_no_stiff =(Wcyl + 2 * Whead) 

    # Size stiffeners
    tanktype = "outer"
    Nstiff = Nmain + Ninterm

    Wmainstiff = stiffener_weight(tanktype, Winnertank / Nmain, Rtank_outer, 
                                s_a, ρouter, Nstiff, θ1, θ2, l_cyl, Eouter) #Weight of one main stiffener
    Wintermstiff = stiffener_weight(tanktype, 0.0, Rtank_outer,
                                s_a, ρouter, Nstiff, θ1, θ2, l_cyl, Eouter) #Weight of one intermediate stiffener

    Wstiff = Nmain * Wmainstiff + Ninterm * Wintermstiff #Total stiffener weight

    Wtank = (Wtank_no_stiff + Wstiff) * (1 + ftankadd)

    return Wtank, Wcyl, Whead, Wstiff, S_outer, Shead, Scyl, t_cyl, t_head
end


