"""
    calculate_thrust_from_ROC!(ac, ip, imission)

A model to estimate the required thrust based on the rate of climb and lift and drag conditions.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object
    - `ip::Int64`: mission point index
    - `imission::Int64`: mission index
     
    **Outputs:**
    No outputs. The `ac` object gets modified with the thrust Fe.
"""
function calculate_thrust_from_ROC!(ac, ip, imission)
    parg, _, para, pare, _, _, _, wing, _, _, _, _ = unpack_ac(ac, imission, ip = ip)

    neng = ac.parg[igneng]
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
end