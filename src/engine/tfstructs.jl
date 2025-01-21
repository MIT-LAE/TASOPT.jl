"""
    flowStation

A type meant to represent the static and stagnation states (and derivatives) of a flow at 
a given internal engine station. Constructor requires `n::Int` as number of non-fuel
air constituents.

# Fields:
- `M::Float64` : Mach number
- `A::Float64` : Area
- `u::Float64` : Speed

- `p?`::Float64 : Pressure; ? = s for static, t for total
- `T?`::Float64 : Temperature; ? = s for static, t for total
- `h?`::Float64 : Enthalpy; ? = s for static, t for total
- `cp?`::Float64 : Specific heat at constant pressure; ? = s for static, t for total
- `cv?`::Float64 : Specific heat at constant volume; ? = s for static, t for total
- `gam?`::Float64 : Specific heat ratio; ? = s for static, t for total
- `a?`::Float64 : Speed of sound; ? = s for static, t for total

- Derivatives of each static and total quantity wrt each other
  
"""
struct flowStation
    # ---- Physical Flow
    M::Float64
    A::Float64
    u::Float64
    u_M::Float64

    # ---- Static Quantities
    ps::Float64
    Ts::Float64
    hs::Float64
    rhos::Float64
    ss::Float64
    cps::Float64
    cvs::Float64
    gams::Float64
    as::Float64

    # ---- Total Quantities
    pt::Float64
    Tt::Float64
    ht::Float64
    rhot::Float64
    st::Float64
    cpt::Float64
    cvt::Float64
    gamt::Float64
    at::Float64

    # ---- Static Derivatives
    ps_ps::Float64
    ps_pt::Float64
    ps_ss::Float64
    ps_st::Float64
    ps_Ts::Float64
    ps_Tt::Float64
    ps_hs::Float64
    ps_ht::Float64
    ps_M::Float64
    ps_pf::Float64
    ps_pl::Float64
    ps_ph::Float64
    ps_mf::Float64
    ps_ml::Float64
    ps_mh::Float64
    ps_ef::Float64
    ps_el::Float64
    ps_eh::Float64

    Ts_ps::Float64
    Ts_pt::Float64
    Ts_ss::Float64
    Ts_st::Float64
    Ts_Ts::Float64
    Ts_Tt::Float64
    Ts_hs::Float64
    Ts_ht::Float64
    Ts_M::Float64
    Ts_pf::Float64
    Ts_pl::Float64
    Ts_ph::Float64
    Ts_mf::Float64
    Ts_ml::Float64
    Ts_mh::Float64
    Ts_ef::Float64
    Ts_el::Float64
    Ts_eh::Float64

    hs_ps::Float64
    hs_pt::Float64
    hs_ss::Float64
    hs_st::Float64
    hs_Ts::Float64
    hs_Tt::Float64
    hs_hs::Float64
    hs_ht::Float64
    hs_M::Float64
    hs_pf::Float64
    hs_pl::Float64
    hs_ph::Float64
    hs_mf::Float64
    hs_ml::Float64
    hs_mh::Float64
    hs_ef::Float64
    hs_el::Float64
    hs_eh::Float64

    rhos_ps::Float64
    rhos_pt::Float64
    rhos_ss::Float64
    rhos_st::Float64
    rhos_Ts::Float64
    rhos_Tt::Float64
    rhos_hs::Float64
    rhos_ht::Float64
    rhos_M::Float64
    rhos_pf::Float64
    rhos_pl::Float64
    rhos_ph::Float64
    rhos_mf::Float64
    rhos_ml::Float64
    rhos_mh::Float64
    rhos_ef::Float64
    rhos_el::Float64
    rhos_eh::Float64

    ss_ps::Float64
    ss_pt::Float64
    ss_ss::Float64
    ss_st::Float64
    ss_Ts::Float64
    ss_Tt::Float64
    ss_hs::Float64
    ss_ht::Float64
    ss_M::Float64
    ss_pf::Float64
    ss_pl::Float64
    ss_ph::Float64
    ss_mf::Float64
    ss_ml::Float64
    ss_mh::Float64
    ss_ef::Float64
    ss_el::Float64
    ss_eh::Float64

    # ---- Total Derivatives
    pt_ps::Float64
    pt_pt::Float64
    pt_ss::Float64
    pt_st::Float64
    pt_Ts::Float64
    pt_Tt::Float64
    pt_hs::Float64
    pt_ht::Float64
    pt_M::Float64
    pt_pf::Float64
    pt_pl::Float64
    pt_ph::Float64
    pt_mf::Float64
    pt_ml::Float64
    pt_mh::Float64
    pt_ef::Float64
    pt_el::Float64
    pt_eh::Float64

    Tt_ps::Float64
    Tt_pt::Float64
    Tt_ss::Float64
    Tt_st::Float64
    Tt_Ts::Float64
    Tt_Tt::Float64
    Tt_hs::Float64
    Tt_ht::Float64
    Tt_M::Float64
    Tt_pf::Float64
    Tt_pl::Float64
    Tt_ph::Float64
    Tt_mf::Float64
    Tt_ml::Float64
    Tt_mh::Float64
    Tt_ef::Float64
    Tt_el::Float64
    Tt_eh::Float64

    ht_ps::Float64
    ht_pt::Float64
    ht_ss::Float64
    ht_st::Float64
    ht_Ts::Float64
    ht_Tt::Float64
    ht_hs::Float64
    ht_ht::Float64
    ht_M::Float64
    ht_pf::Float64
    ht_pl::Float64
    ht_ph::Float64
    ht_mf::Float64
    ht_ml::Float64
    ht_mh::Float64
    ht_ef::Float64
    ht_el::Float64
    ht_eh::Float64

    rhot_ps::Float64
    rhot_pt::Float64
    rhot_ss::Float64
    rhot_st::Float64
    rhot_Ts::Float64
    rhot_Tt::Float64
    rhot_hs::Float64
    rhot_ht::Float64
    rhot_M::Float64
    rhot_pf::Float64
    rhot_pl::Float64
    rhot_ph::Float64
    rhot_mf::Float64
    rhot_ml::Float64
    rhot_mh::Float64
    rhot_ef::Float64
    rhot_el::Float64
    rhot_eh::Float64

    st_ps::Float64
    st_pt::Float64
    st_ss::Float64
    st_st::Float64
    st_Ts::Float64
    st_Tt::Float64
    st_hs::Float64
    st_ht::Float64
    st_M::Float64
    st_pf::Float64
    st_pl::Float64
    st_ph::Float64
    st_mf::Float64
    st_ml::Float64
    st_mh::Float64
    st_ef::Float64
    st_el::Float64
    st_eh::Float64


    # Constructor with `n` parameter for `_al` vector sizes
    function flowStation( 
        M::Float64 = 0.0, A::Float64 = 1.0, u::Float64 = 0.0,
        ps::Float64 = 1.0, Ts::Float64 = 300.0, hs::Float64 = 1.0, rhos::Float64 = 0.0, ss::Float64 = 0.0,
        cps::Float64 = 1.0, cvs::Float64 = 1.0, gams::Float64 = 1.4, as::Float64 = 340.0,
        pt::Float64 = 1.0, Tt::Float64 = 300.0, ht::Float64 = 1.0, rhot::Float64 = 0.0, st::Float64 = 0.0,
        cpt::Float64 = 1.0, cvt::Float64 = 1.0, gamt::Float64 = 1.4, at::Float64 = 340.0
    )
        return new(
            # ---- Physical Flow
            M, A, u, 0.0,
    
            # ---- Static Quantities
            ps, Ts, hs, rhos, ss, cps, cvs, gams, as,
    
            # ---- Total Quantities
            pt, Tt, ht, rhot, st, cpt, cvt, gamt, at,
    
            # ---- Static Derivatives
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    
            # ---- Total Derivatives
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
    end


    function flowStation(
        M::Float64 = 0.0, A::Float64 = 1.0, u::Float64 = 0.0,
        ps::Float64 = 1.0, Ts::Float64 = 300.0, hs::Float64 = 1.0, ss::Float64 = 0.0,
        cps::Float64 = 1.0, cvs::Float64 = 1.0, gams::Float64 = 1.4, as::Float64 = 340.0,
        pt::Float64 = 1.0, Tt::Float64 = 300.0, ht::Float64 = 1.0, st::Float64 = 0.0,
        cpt::Float64 = 1.0, cvt::Float64 = 1.0, gamt::Float64 = 1.4, at::Float64 = 340.0
    )
        return new(
            # ---- Physical Flow
            M, A, u,

            # ---- Static Quantities
            ps, Ts, hs, ss, cps, cvs, gams, as,

            # ---- Total Quantities
            pt, Tt, ht, st, cpt, cvt, gamt, at,
            
            # ---- Static Derivatives
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            
            # ---- Total Derivatives
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
    end
end # flowStation


"""
    compressor

A type meant to store information relevant to an axial compressor. Constructor
requires design quantities for the device.

# Fields:
- `prD::Float64` : Design pressure ratio
- `mbD::Float64` : Design corrected mass flow
- `epD::Float64` : Design polytropic efficiency
- `NbD::Float64` : Design corrected speed

- `map::compressorTbl` : compressorTbl representation of the associated map
- `oob_map::SVector{9, Float64}` : static array of size 9 containing drela map coefficients
                                   used for out-of-bounds map calcualtions

- `pr::Float64` : pressure ratio
- `mb::Float64` : corrected mass flow
- `ep::Float64` : polytropic efficiency
- `Nb::Float64` : corrected speed

- `Nb_?::Float64` : Derivative corrected speed
- `ep_?::Float64` : Derivative polytropic efficiency

"""
struct compressor
# ---- Design values
prD::Float64
mbD::Float64
epD::Float64
NbD::Float64

# ---- Map (and out-of-bounds coefficients)
map::compressorTbl
oob_map::SVector{9, Float64}

# ---- Current values
pr::Float64
mb::Float64
ep::Float64
Nb::Float64

# ---- Derivatives
Nb_pr::Float64
Nb_mb::Float64

ep_pr::Float64
ep_mb::Float64

function compressor(prD::Float64, mbD::Float64, epD::Float64, NbD::Float64,
    map::compressorTbl, oob_map::SVector{9, Float64}; pr::Float64=0., 
    mb::Float64=0., ep::Float64=0., Nb::Float64=0.
)

    if pr == 0.
        pr = prD
    end
    if mb == 0.
        mb = mbD
    end
    if ep == 0.
        ep = epD
    end
    if Nb == 0.
        Nb = NbD
    end

    Nb_pr = 0.0
    Nb_mb = 0.0
    ep_pr = 0.0
    ep_mb = 0.0

    new(prD, mbD, epD, NbD, map, oob_map, pr, mb, ep, Nb, Nb_pr, Nb_mb, ep_pr, ep_mb)
end

end # compressor


"""
    turbine

A type meant to store information relevant to an axial compressor. Constructor
requires design quantities for the device.

# Fields:
- `prD::Float64` : Design pressure ratio
- `mbD::Float64` : Design corrected mass flow
- `epD::Float64` : Design polytropic efficiency
- `NbD::Float64` : Design corrected speed

- `map::SVector{2, Float64}` : Drela turbine coefficients

- `pr::Float64` : pressure ratio
- `mb::Float64` : corrected mass flow
- `ep::Float64` : polytropic efficiency
- `Nb::Float64` : corrected speed

- `ep_?::Float64` : Derivative polytropic efficiency

"""
struct turbine
# ---- Design values
prD::Float64
mbD::Float64
epD::Float64
NbD::Float64

# ---- Map (and out-of-bounds coefficients)
map::SVector{2, Float64}

# ---- Current values
pr::Float64
mb::Float64
ep::Float64
Nb::Float64

# ---- Derivatives
ep_dh::Float64
ep_mb::Float64
ep_Nb::Float64
ep_Tt41::Float64
ep_cpt41::Float64
ep_Rt41::Float64

function compressor(prD::Float64, mbD::Float64, epD::Float64, NbD::Float64,
    map::SVector{2, Float64}; pr::Float64=0., mb::Float64=0., ep::Float64=0.,
    Nb::Float64=0.
)

    if pr == 0.
        pr = prD
    end
    if mb == 0.
        mb = mbD
    end
    if ep == 0.
        ep = epD
    end
    if Nb == 0.
        Nb = NbD
    end

    ep_dh = 0.0
    ep_mb = 0.0
    ep_Nb = 0.0
    ep_Tt41 = 0.0
    ep_cpt41 = 0.0
    ep_Rt41 = 0.0

    new(prD, mbD, epD, NbD, map, pr, mb, ep, Nb, ep_dh, ep_mb, ep_Nb, ep_Tt41, ep_cpt41, ep_Rt41)
end

end # turbine