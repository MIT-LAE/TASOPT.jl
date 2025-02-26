function constant_TSFC_engine!(ac, case::String, imission::Int64, ip::Int64, initializes_engine::Bool, iterw::Int64 = 0)
    #Unpack data storage arrays
    pari, parg, _, para, pare, _, _, wing, _, _, _ = unpack_ac(ac, imission, ip = ip)
    TSFC = pare[ieTSFC]
    neng = parg[igneng]

    if (ip in range(ipstatic, ipclimbn)) #If aircraft is in climb, thrust is determined from the climb rate
        
        WMTO = parg[igWMTO]
        S = wing.layout.S
        DoL = para[iaCD] / para[iaCL]
        W = para[iafracW] * WMTO
        BW = W + para[iaWbuoy]
        CL = para[iaCL]
        ρ = pare[ierho0]

        ROC = para[iaROCdes]

        #Solve with numerical & analytical solution
        f(γ) = ROC - sin(γ) * sqrt(2*BW*cos(γ) / (ρ*S*CL))
        γ = find_zero(f, para[iagamV]) #Solve for climb angle
        A = sin(γ)
        B = DoL
        ϕ = (sqrt(-A^2*B^6 - 2*A^2*B^4 - A^2*B^2 + B^6 + 2*B^4 + B^2) + A*B^2 + A)/(B^2 + 1)
        Ftotal = BW * ϕ
        Fe = Ftotal / neng #required thrust per engine
        pare[ieFe] = Fe

    else
        Fe = pare[ieFe]
        
    end

    mfuel_per_eng =  TSFC * Fe / gee #Calculate fuel flow per engine from set TSFC
    pare[iemfuel] = mfuel_per_eng * neng
end