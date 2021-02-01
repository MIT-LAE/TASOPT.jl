# Properties of hydrogen at saturation
# Saturated liquid
    ρl(P) = 0.3175*P^2 - 4.1121*P + 74.544 # Fit to above data P is in atm!
    ul(P) = -2.8443*P^2 + 34.84P  - 33.044
    hl(P) = -2.7899*P^2 + 36.244P - 33.072
# Saturated vapor
    ρg(P) =  9.0e-5*P^2 + 1.1838P + 0.1589
    ug(P) = -2.3398*P^2 + 13.435P + 359.02
    hg(P) = -3.7721*P^2 + 22.094P + 427.65


# Let yg be the ullage volume fraction - 3% ullage ⟹ a volume fill factor of 97%
ρmix(yg,P) = ρl(P)*(1-yg) + ρg(P)*yg

quality(yg, P) = 1/( 1 + (ρl(P)/ρg(P)) * (1-yg)/yg )

function φ(P,δP,yg)
    ρ = 0.5*(ρmix(yg, P) + ρmix(yg, P+δP))
    x = ( 1/ρ - 1/ρl(P) ) / ( 1/ρg(P) - 1/ρl(P) )

    uᵢ   = x*ug(P)    + (1-x)*ul(P)
    uᵢ₊₁ = x*ug(P+δP) + (1-x)*ul(P+δP)
    φ = (ρ*(uᵢ₊₁ - uᵢ)/(δP*101.325))^-1

    return φ
end

function dPdt_LH2(ΔT, mdot_out, P, yg)
    Req = 0.225
    dqdt = 1.3*ΔT/Req
    Vfuel = 200.0
    ρbar = ρg(P)/(ρg(P) - ρl(P))

    dPdt = φ(P, 0.01, yg)/Vfuel * (dqdt - mdot_out*(hl(P)-hg(P))*(ρbar))

    return dPdt
end

mdot = [0.617606, 0.52916, 0.437406, 0.359918, 0.295704, 0.291292, 0.271806]
time = [174.10271765493297, 340.5402718139433, 584.5015556778529, 1155.4673259995732, 1155.4673259995732, 22678.920538191625, 22678.920538191625]

P = zero(mdot)
yg = zero(mdot)
P[1] = 1.2
yg[1] = 0.1
ΔT = 226.0 - 20.0
for (i,t) in enumerate(time)
    
    if(i<7)
        yg[i+1] = yg[i] + mdot[i]/ρmix(yg[i], P[i])
        
        P[i+1] = P[i] + dPdt_LH2(ΔT, mdot[i], P[i], yg[i])*(time[i+1] - time[i])
        println(dPdt_LH2(ΔT, mdot[i], P[i], yg[i]))
    end
end
P