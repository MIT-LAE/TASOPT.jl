""" 
    constant_TSFC_engine!(ac, d1, imission::Int64, ip::Int64, d2, d3)

A very simple model of an aircraft engine consuming some fuel. The fuel flow is determined by the thrust required 
by the aircraft. The thrust is calculated based on the climb rate of the aircraft. The fuel flow is calculated using 
the thrust specific fuel consumption (TSFC) of the engine, which is an input.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object
    - `imission::Int64`: mission index
    - `ip::Int64`: mission point index
    - `d`s: dummy variables used for compatibility with other sizing functions

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine parameters.
"""
function constant_TSFC_engine!(ac, d1, imission::Int64, ip::Int64, d2, d3 = 0)
    #Unpack data storage arrays
    pari, parg, _, para, pare, _, _, wing, _, _, _ = unpack_ac(ac, imission, ip = ip)
    TSFC = pare[ieTSFC] #Extract TSFC at this mission point
    neng = parg[igneng]

    if (ip in range(ipstatic, ipclimbn)) #If aircraft is in climb, thrust is determined from the climb rate
        
        #Calculate climb angle from desired climb rate
        WMTO = parg[igWMTO]
        S = wing.layout.S
        DoL = para[iaCD] / para[iaCL]
        W = para[iafracW] * WMTO
        BW = W + para[iaWbuoy]
        CL = para[iaCL]
        ρ = pare[ierho0]

        ROC = para[iaROCdes] #desired rate of climb in m/s

        #Solve with numerical & analytical solution
        f(γ) = ROC - sin(γ) * sqrt(2*BW*cos(γ) / (ρ*S*CL))
        γ = find_zero(f, para[iagamV]) #Solve for climb angle
        A = sin(γ)
        B = DoL
        ϕ = (sqrt(-A^2*B^6 - 2*A^2*B^4 - A^2*B^2 + B^6 + 2*B^4 + B^2) + A*B^2 + A)/(B^2 + 1)#Closed form solution
        Ftotal = BW * ϕ #Total thrust required for climb
        Fe = Ftotal / neng #required thrust per engine
        pare[ieFe] = Fe #Store computed thrust

    else
        Fe = pare[ieFe] #In other mission points, thrust is already computed
        
    end

    mfuel_per_eng =  TSFC * Fe / gee #Calculate fuel flow per engine from set TSFC
    pare[iemfuel] = mfuel_per_eng * neng #Store total fuel flow for all engines
end