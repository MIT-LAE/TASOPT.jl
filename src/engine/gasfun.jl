using Zygote

include("gasdata.jl")

"""
    gasfun(igas, t)

Computes properties of a thermally-perfect gas
with some variable specific heat cp[T].

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `igas`: index specifying the gas (see if blocks below for list)
    - `t`: temperature T in Kelvin

    **Output:**
    - `s`: entropy-complement function s[T]
    - `s_t`: ds/dT
    - `h`: complete enthalpy function h[T]
    - `h_t`: dh/dT
    - `cp`: specific heat cp[T]
    - `r`: ideal-gas constant R


The adiabatic pressure change over a process 1->2 with some polytropic efficiency epol is:

``\\ p2  = \\ p1  exp [   epol   (s2-s1)/R ] ``   compression

``\\ p2  = \\ p1  exp [ (1/epol) (s2-s1)/R ] ``   expansion
"""
function gasfun(igas, t)

    if (igas == 1)
        s, h, cp, r = gas_N2(t)
    elseif (igas == 2)
        s, h, cp, r = gas_O2(t)
    elseif (igas == 3)
        s, h, cp, r = gas_CO2(t)
    elseif (igas == 4)
        s, h, cp, r = gas_H2O(t)
    elseif (igas == 5)
        s, h, cp, r = gas_Ar(t)
    elseif (igas == 11)
        s, h, cp, r = gas_CH4(t)
    elseif (igas == 12)
        s, h, cp, r = gas_C2H6(t)
    elseif (igas == 13)
        s, h, cp, r = gas_C3H8(t)
    elseif (igas == 14)
        s, h, cp, r = gas_C4H10(t)
    elseif (igas == 18)
        s, h, cp, r = gas_C8H18(t)
    elseif (igas == 24)
        s, h, cp, r = gas_C14H30(t)
    elseif (igas == 40)
        s, h, cp, r = gas_H2(t)
    else
        error("GASFUN: undefined gas index: ", igas, " & code=0")
    end
    s_t = cp / t
    h_t = cp
    return s, s_t, h, h_t, cp, r
end # gasfun

"""
    gaschem(igas)

Returns number of C,H,O,N atoms in gas molecule, for the gases implemented in function gasfun above.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `igas`: index specifying the gas (see if blocks below for list)

    **Output:**
    - `nchon(.)`: number of C,H,O,N atoms in gas molecule

"""
function gaschem(igas)

    kc = 1
    kh = 2
    ko = 3
    kn = 4

    nchon = [0, 0, 0, 0]
    buf = Zygote.Buffer(nchon, length(nchon))
    for i = 1:length(nchon)
        buf[i] = nchon[i]
    end

    if (igas == 1)  #   N2
        buf[kn] = 2
    elseif (igas == 2)  #   O2
        buf[ko] = 2
    elseif (igas == 3)  #   CO2
        buf[kc] = 1
        buf[ko] = 2
    elseif (igas == 4)  #   H2O
        buf[kh] = 2
        buf[ko] = 1
    elseif (igas == 5)  #   Ar

    elseif (igas == 11)  #   CH4
        buf[kc] = 1
        buf[kh] = 4
    elseif (igas == 12)  #   C2H6
        buf[kc] = 2
        buf[kh] = 6
    elseif (igas == 13)  #   C3H8
        buf[kc] = 3
        buf[kh] = 8
    elseif (igas == 14)  #   C4H10
        buf[kc] = 4
        buf[kh] = 10
    elseif (igas == 18)  #   C8H18
        buf[kc] = 8
        buf[kh] = 18
    elseif (igas == 24)  #   C14H30
        buf[kc] = 14
        buf[kh] = 30
    elseif (igas == 40)  #   H2
        buf[kc] = 0
        buf[kh] = 2
    else
        error("GASFUN: undefined gas index: ", igas, " & code=0")
    end
    nchon = copy(buf)
    return nchon
end # gaschem

function gas_N2(t1, t, tl, cp, cpt, h, s)
    r = 296.94
    hform = 0.0000
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_N2

function gas_O2(t1, t, tl, cp, cpt, h, s)
    r = 259.82
    hform = 0.000
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_O2

function gas_Ar(t1, t, tl, cp, cpt, h, s)
    r = 208.00
    hform = 0.0000
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_Ar

function gas_CO2(t1, t, tl, cp, cpt, h, s)
    r = 188.96
    hform = -0.89430E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_CO2

function gas_H2O(t1, t, tl, cp, cpt, h, s)
    r = 461.91
    hform = -0.13433E+08
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_H2O

function gas_CH4(t1, t, tl, cp, cpt, h, s)
    r = 519.65
    hform = -0.46750E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_CH4

function gas_C2H6(t1, t, tl, cp, cpt, h, s)
    r = 277.15
    hform = -0.27930E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_C2H6

function gas_C3H8(t1, t, tl, cp, cpt, h, s)
    r = 188.50
    hform = -0.23570E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_C3H8

function gas_C4H10(t1, t, tl, cp, cpt, h, s)
    r = 143.30
    hform = -0.21640E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_C4H10

function gas_C8H18(t1, t, tl, cp, cpt, h, s)
    r = 72.900
    hform = -0.18280E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_C8H18

function gas_C14H30(t1, t, tl, cp, cpt, h, s)
    r = 167.00
    hform = -0.16700E+07
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
    # gas_C14H30 
end

function gas_H2(t1, t, tl, cp, cpt, h, s)
    r = 4124.9
    hform = 0.0000
    get_thermo(r, hform, t1, t, tl, cp, cpt, h, s)
end # gas_H2

macro define_gas_method(func_name, gas_name)
    func_symbol = esc(func_name)
    gas_symbol = esc(gas_name)
    return quote
        $func_symbol(t1) = $func_symbol(
            t1,
            $gas_symbol.t,
            $gas_symbol.tl,
            $gas_symbol.cp,
            $gas_symbol.cpt,
            $gas_symbol.h,
            $gas_symbol.s,
        )
    end
end

@define_gas_method gas_N2 N2
@define_gas_method gas_O2 O2
@define_gas_method gas_Ar Ar
@define_gas_method gas_CO2 CO2
@define_gas_method gas_H2O H2O
@define_gas_method gas_CH4 CH4
@define_gas_method gas_C2H6 C2H6
@define_gas_method gas_C3H8 C3H8
@define_gas_method gas_C4H10 C4H10
@define_gas_method gas_C8H18 C8H18
@define_gas_method gas_C14H30 C14H30
@define_gas_method gas_H2 H2

"""
    findsegment(x::T, xarr::Vector{T})

Uses bisection to find the right interval of the array `xarr` where x lies.
Returns im and io s.t. `xarr[im] < x < xarr[io]`.

Additionally returns the interval `dx = xarr[io] - xarr[im]`
"""
function findsegment(x::Float64, xarr::AbstractArray)

    if isnan(x)
        error("Oops, you're searching for a NaN! Go fix your bug!")
    end

    io::Int = length(xarr)

    if x ≤ xarr[1]
        im = 1
        io = 2
    elseif x ≥ xarr[end]
        im = io - 1
    else
        im = searchsortedlast(xarr, x)
        io = im + 1
    end

    dx = xarr[io] - xarr[im]
    t = (x - xarr[im]) / dx

    return io, im, dx, t
end

"""
Helper function that returns the thermo properties by table interpolation
"""
function get_thermo(
    r::T,
    hform::T,
    t1::T,
    t::Vector{T},
    tl::Vector{T},
    cp::Vector{T},
    cpt::Vector{T},
    h::Vector{T},
    s::Vector{T},
) where {T<:AbstractFloat}

    i, im, dt, f = findsegment(t1, t)

    dtl = tl[i] - tl[im]
    fl = (log(t1) - tl[im]) / dtl
    A = 1.0 - f
    B = 1.0 - fl
    #- - - - - - - - - - -
    s1 =
        B * s[im] +
        fl * s[i] +
        fl * B * (B * (dtl * cp[im] - s[i] + s[im]) - fl * (dtl * cp[i] - s[i] + s[im]))

    h1 =
        hform +
        A * h[im] +
        f * h[i] +
        f * A * (A * (dt * cp[im] - h[i] + h[im]) - f * (dt * cp[i] - h[i] + h[im]))

    cp1 =
        A * cp[im] +
        f * cp[i] +
        f * A * (A * (dt * cpt[im] - cp[i] + cp[im]) - f * (dt * cpt[i] - cp[i] + cp[im]))

    r1 = r
    #- - - - - - - - - - -
    return s1, h1, cp1, r1
end  # function get_thermo
