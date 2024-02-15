function freestream_heat_coeff(z, M, xftank, Tw)
    #Use ISA function to calculate freestream conditions
    Tair, p, _, a, _ = atmos(z / 1e3)
    u = M * a #freestrean velocity

    #Assumed parameters for air
    Pr = 0.71
    γ = 1.4
    cp = 1005
    R = 287
    
    #Find h for air
    # This code uses the reference temperature model and the Chilton-Colburn analogy
    T_s = Tair * (0.5 * (1 + Tw/Tair) + 0.16 * Pr^(1/3) * (γ - 1)/2 * M^2) #Reference temperature
    
    #Find viscosity from Sutherland's law
    μ0 = 1.716e-5
    S_μ = 111
    T0 = 273
    μ_s = μ0 * (T_s / T0)^(3 / 2) * ( (T0 + S_μ) / (T_s + S_μ) )

    ρ_s = p / (R * T_s) #density at reference temperature

    Re_xftank = ρ_s * u * xftank / μ_s
    cf_xftank = 0.02296 / (Re_xftank)^0.139
    #Calculate Stanton number using Reynolds analogy
    St_air = cf_xftank / (2 * Pr^(2/3)) #Chilton-Colburn analogy
    hconvair = St_air * ρ_s *u* cp #In W/(m^2 K)

    return hconvair, Tair
end
