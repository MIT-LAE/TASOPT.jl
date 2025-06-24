"""
$TYPEDEF

Simplified model of a stranded (i.e., Litz) wire  with an outer insulation layer.
Default values are for copper based on https://pubs.aip.org/aip/jpr/article/8/4/1147/242248/Electrical-resistivity-of-copper-gold-palladium 

# Fields
$TYPEDFIELDS

"""
@kwdef mutable struct Cable
    """Conductor"""
    cond::Conductor = Conductor("Cu")
    """Insulator"""
    ins::Insulator = Insulator("polyamide")
    """Temperature of conductors [K]"""
    Tcon::Float64 = 273.1 + 50.0 #defaults to 50ᵒC
    """Max current density [A/m²]"""
    Jmax::Float64 = 5e6 # 5 A/mm² = 5e6 A/m²
    """Packing factor [-]"""
    kpf::Float64 = 0.91
    """Length [m]"""
    l::Float64 = NaN
    """Mass [kg]"""
    mass::Float64 = NaN
    """Weight [N]"""
    W::Float64 = NaN
    """Resistance [Ω]"""
    R::Float64 = NaN
    """Voltage [V]"""
    V::Float64 = NaN
end

"""
    (c::Cable)(P::Float64, V::Float64, lcable::Float64)

Sizes the cable given a power, voltage, and length of the cable and returns 
the *efficiency function* for the sized cable.
"""
function (c::Cable)(P::Float64, V::Float64, lcable::Float64)
    cond = c.cond
    ins  = c.ins
    c.V = V
    c.l = lcable
    # Calculate current
    I = P/V
    
    # Get resistivity at given temperature 
    ρcon = resistivity(cond, c.Tcon)
    
    # Geometry
    tins = V/ins.Emax
    Acon = I/c.Jmax
    ri = sqrt(Acon/(π*c.kpf)) # A = πr²×kpf
    ro = ri + tins
    Ains = π*(ro^2 - ri^2)
    c.mass = lcable*(Acon*cond.ρ + Ains*ins.ρ)
    c.W = c.mass*gee
    
    # Resistance
    c.R = ρcon * lcable/Acon

    # Loss function
    """
    Returns the efficiency of the cable for the given power input assuming constant voltage.
    """
    function f(Pin::Float64)
        
        Iin = Pin/c.V

        if Pin>P
            @warn """Power ($(Pin/1e3) kW) exceeds design power for cable ($(P/1e3) kW).
             Current density = $(Iin/Acon) > Jmax ($(c.Jmax)) A/m²"""
        end

        Lohmic = Iin^2 * c.R # = I²R loss
        Pout = Pin - Lohmic
        η = Pout/Pin
        return η
    end  # function 
    return f
end