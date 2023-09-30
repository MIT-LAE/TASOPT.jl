module atmosphere

export atmos

"""
    atmos(h)
    
Atmospheric functions ` T(h)`, `ρ(h)` etc
valid to `h`=20km, `p(h)` valid to `h`=70km.

Also calculates viscosity using Sutherland's law.

Units:
- [h]   = km ASL
- [T]   = Kelvin
- [p]   = Pa
- [ρ]   = kg/m^3
- [a]   = m/s
- [μ]   = kg/m-s 
"""
function atmos(h)

 pSL = 1.0132e5 # Pa
 TSL    = 288.2    # K
 Tpause = 216.65   # K
 Tblend = 2.0      # K (tropopause blending T range)
 Tlapse = -6.5     # K/km
 
 μSL  = 1.78e-5  # kg/m-s
 Tsuth = 110.0    # K
 
 cp  = 1004.0   # J/kg-K
 ɣ = 1.4 
 
 T = log(1.0+exp((TSL+Tlapse*h-Tpause)/Tblend))*Tblend + Tpause
 p = pSL*exp( - 0.11800*h /(1.0 + 0.0020*h) - 0.00198*h^2/(1.0 + 0.0006*h^2) )
 ρ = ɣ*p/((ɣ-1.0)*cp*T)
 a = sqrt(ɣ*p/ρ)
 μ = μSL * sqrt(T/TSL)^3 * (TSL+Tsuth)/(T+Tsuth)
 
 return T,p,ρ,a,μ

end # atmos

end