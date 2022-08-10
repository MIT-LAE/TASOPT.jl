# Properties of hydrogen at saturation from NIST
# NIST data for saturated LH2:
# https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=B5000001&Type=SatT&Digits=5&PLow=1&PHigh=3&PInc=0.05&RefState=DEF&TUnit=K&PUnit=atm&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm
    # P is in [atm] 
    # Saturated liquid
        # Density
        ρl(P)  =   0.3175*P^2 - 4.1121*P + 74.544 # returns in kg/m3
        ρl′(P) = (2*0.3175*P   - 4.1121)/101325   # returns in kg/m3/Pa
        
        # Internal energy
        # ul(P) = -2.8443*P^2 + 34.84P  - 33.044 # returns in kJ/kg
        ul(P)  =   0.8443*P^3 -   7.91*P^2 + 44.44*P - 38.735
        ul′(P) = 0.8443*3*P^2 - 2*7.91*P   + 44.44            # returns kJ/kg/atm

        hl(P) = -2.7899*P^2 + 36.244P - 33.072 # returns in kJ/kg
        
    # Saturated vapor
        # Density
        ρg(P)  =    9.0e-5*P^2 + 1.1838P + 0.1589 # returns in kg/m3
        ρg′(P) =  (2*9.0e-5*P   + 1.1838)/101325           # returns in kg/m3/Pa
        
        # Internal Energy
        # ug(P) = -2.3398*P^2 + 13.435P + 359.02 # returns in kJ/kg
        ug(P)  =   0.674*P^3 -   6.3839*P^2 + 21.099*P + 354.47    # returns in kJ/kg
        ug′(P) = 3*0.674*P^2 - 2*6.3839*P   + 21.099    # returns in kJ/kg/atm
        
        hg(P) = -3.7721*P^2 + 22.094P + 427.65 # returns in kJ/kg

    # mixed
        umix(x,P) = x*ug(P) + (1-x)*ul(P)
        # want in kJ/kg / kPa i.e. unitless so ÷ by 101.325
        # u′mix(x,P) = (x*(-2.3398*2*P + 13.435) +(1-x)*(-2.8443*P*2 + 34.84)) /101.325 
        u′mix(x,P) = (x*ug′(P) + (1-x)*ul′(P) )/101.325 


# Let yg be the ullage volume fraction - 3% ullage ⟹ a volume fill factor of 97%
ρmix(yg,P) = ρl(P)*(1-yg) + ρg(P)*yg
ρmix_qual(x,P) = (x/ρg(P) + (1-x)/ρl(P))^-1

quality(yg, P) = 1/( 1 + (ρl(P)/ρg(P))*((1-yg)/yg) )

#P in atm returns φ in kJ/kg/kPa = [-]
function φ(P,yg)
    # x = quality(yg, P)
    ρ = ρmix(yg,P) 
    x = ( 1/ρ - 1/ρl(P) ) / ( 1/ρg(P) - 1/ρl(P) )
    # δP = 0.01
    # uᵢ   = x*ug(P)    + (1-x)*ul(P)
    # uᵢ₊₁ = x*ug(P+δP) + (1-x)*ul(P+δP)
    # ∂u_∂P = (uᵢ₊₁ - uᵢ)/(δP*101.325) 
    # φ = (ρ*∂u_∂P)^-1
    φ = (ρ*u′mix(x,P))^-1

    return φ
end

function dPdt_LH2(ΔT, mdot_out, P, yg)
    Req = 0.225
    dqdt = 1.3*ΔT/Req
    Vfuel = 200.0
    ρbar = ρg(P)/(ρg(P) - ρl(P))

    dPdt = φ(P, yg)/Vfuel * (dqdt - mdot_out*(hl(P)-hg(P))*1000*(ρbar))

    return dPdt/101325
end

# Self pressurization 
dPdt_sp(Q, V, P, yg) = φ(P, yg)/V * Q
