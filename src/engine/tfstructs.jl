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
    Mi::Float64
    A::Float64
    u::Float64
    u_Mi::Float64

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
    ps_Mi::Float64
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
    Ts_Mi::Float64
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
    hs_Mi::Float64
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
    rhos_Mi::Float64
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
    ss_Mi::Float64
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
    pt_Mi::Float64
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
    Tt_Mi::Float64
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
    ht_Mi::Float64
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
    rhot_Mi::Float64
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
    st_Mi::Float64
    st_pf::Float64
    st_pl::Float64
    st_ph::Float64
    st_mf::Float64
    st_ml::Float64
    st_mh::Float64
    st_ef::Float64
    st_el::Float64
    st_eh::Float64

    function flowStation(;
        Mi::Float64 = 0.0, A::Float64 = 1.0, u::Float64 = 0.0, u_Mi::Float64 = 0.0,
        ps::Float64 = 1.0, Ts::Float64 = 0.0, hs::Float64 = 1.0, rhos::Float64 = 0.0, ss::Float64 = 0.0,
        cps::Float64 = 1.0, cvs::Float64 = 1.0, gams::Float64 = 1.0, as::Float64 = 0.0,
        pt::Float64 = 1.0, Tt::Float64 = 0.0, ht::Float64 = 1.0, rhot::Float64 = 0.0, st::Float64 = 0.0,
        cpt::Float64 = 1.0, cvt::Float64 = 1.0, gamt::Float64 = 1.0, at::Float64 = 0.0,

        ps_ps::Float64 = 0.0, ps_pt::Float64 = 0.0, ps_ss::Float64 = 0.0, ps_st::Float64 = 0.0, 
        ps_Ts::Float64 = 0.0, ps_Tt::Float64 = 0.0, ps_hs::Float64 = 0.0, ps_ht::Float64 = 0.0, 
        ps_Mi::Float64 = 0.0, ps_pf::Float64 = 0.0, ps_pl::Float64 = 0.0, ps_ph::Float64 = 0.0, 
        ps_mf::Float64 = 0.0, ps_ml::Float64 = 0.0, ps_mh::Float64 = 0.0, ps_ef::Float64 = 0.0, 
        ps_el::Float64 = 0.0, ps_eh::Float64 = 0.0,
        
        Ts_ps::Float64 = 0.0, Ts_pt::Float64 = 0.0, Ts_ss::Float64 = 0.0, Ts_st::Float64 = 0.0, 
        Ts_Ts::Float64 = 0.0, Ts_Tt::Float64 = 0.0, Ts_hs::Float64 = 0.0, Ts_ht::Float64 = 0.0, 
        Ts_Mi::Float64 = 0.0, Ts_pf::Float64 = 0.0, Ts_pl::Float64 = 0.0, Ts_ph::Float64 = 0.0, 
        Ts_mf::Float64 = 0.0, Ts_ml::Float64 = 0.0, Ts_mh::Float64 = 0.0, Ts_ef::Float64 = 0.0, 
        Ts_el::Float64 = 0.0, Ts_eh::Float64 = 0.0,

        hs_ps::Float64 = 0.0, hs_pt::Float64 = 0.0, hs_ss::Float64 = 0.0, hs_st::Float64 = 0.0, 
        hs_Ts::Float64 = 0.0, hs_Tt::Float64 = 0.0, hs_hs::Float64 = 0.0, hs_ht::Float64 = 0.0, 
        hs_Mi::Float64 = 0.0, hs_pf::Float64 = 0.0, hs_pl::Float64 = 0.0, hs_ph::Float64 = 0.0, 
        hs_mf::Float64 = 0.0, hs_ml::Float64 = 0.0, hs_mh::Float64 = 0.0, hs_ef::Float64 = 0.0, 
        hs_el::Float64 = 0.0, hs_eh::Float64 = 0.0, 
        
        rhos_ps::Float64 = 0.0, rhos_pt::Float64 = 0.0, rhos_ss::Float64 = 0.0, rhos_st::Float64 = 0.0, 
        rhos_Ts::Float64 = 0.0, rhos_Tt::Float64 = 0.0, rhos_hs::Float64 = 0.0, rhos_ht::Float64 = 0.0, 
        rhos_Mi::Float64 = 0.0, rhos_pf::Float64 = 0.0, rhos_pl::Float64 = 0.0, rhos_ph::Float64 = 0.0, 
        rhos_mf::Float64 = 0.0, rhos_ml::Float64 = 0.0, rhos_mh::Float64 = 0.0, rhos_ef::Float64 = 0.0, 
        rhos_el::Float64 = 0.0, rhos_eh::Float64 = 0.0, 
        
        ss_ps::Float64 = 0.0, ss_pt::Float64 = 0.0, ss_ss::Float64 = 0.0, ss_st::Float64 = 0.0, 
        ss_Ts::Float64 = 0.0, ss_Tt::Float64 = 0.0, ss_hs::Float64 = 0.0, ss_ht::Float64 = 0.0, 
        ss_Mi::Float64 = 0.0, ss_pf::Float64 = 0.0, ss_pl::Float64 = 0.0, ss_ph::Float64 = 0.0, 
        ss_mf::Float64 = 0.0, ss_ml::Float64 = 0.0, ss_mh::Float64 = 0.0, ss_ef::Float64 = 0.0, 
        ss_el::Float64 = 0.0, ss_eh::Float64 = 0.0, 

        pt_ps::Float64 = 0.0, pt_pt::Float64 = 0.0, pt_ss::Float64 = 0.0, pt_st::Float64 = 0.0, 
        pt_Ts::Float64 = 0.0, pt_Tt::Float64 = 0.0, pt_hs::Float64 = 0.0, pt_ht::Float64 = 0.0, 
        pt_Mi::Float64 = 0.0, pt_pf::Float64 = 0.0, pt_pl::Float64 = 0.0, pt_ph::Float64 = 0.0, 
        pt_mf::Float64 = 0.0, pt_ml::Float64 = 0.0, pt_mh::Float64 = 0.0, pt_ef::Float64 = 0.0, 
        pt_el::Float64 = 0.0, pt_eh::Float64 = 0.0,
        
        Tt_ps::Float64 = 0.0, Tt_pt::Float64 = 0.0, Tt_ss::Float64 = 0.0, Tt_st::Float64 = 0.0, 
        Tt_Ts::Float64 = 0.0, Tt_Tt::Float64 = 0.0, Tt_hs::Float64 = 0.0, Tt_ht::Float64 = 0.0, 
        Tt_Mi::Float64 = 0.0, Tt_pf::Float64 = 0.0, Tt_pl::Float64 = 0.0, Tt_ph::Float64 = 0.0, 
        Tt_mf::Float64 = 0.0, Tt_ml::Float64 = 0.0, Tt_mh::Float64 = 0.0, Tt_ef::Float64 = 0.0, 
        Tt_el::Float64 = 0.0, Tt_eh::Float64 = 0.0,

        ht_ps::Float64 = 0.0, ht_pt::Float64 = 0.0, ht_ss::Float64 = 0.0, ht_st::Float64 = 0.0, 
        ht_Ts::Float64 = 0.0, ht_Tt::Float64 = 0.0, ht_hs::Float64 = 0.0, ht_ht::Float64 = 0.0, 
        ht_Mi::Float64 = 0.0, ht_pf::Float64 = 0.0, ht_pl::Float64 = 0.0, ht_ph::Float64 = 0.0, 
        ht_mf::Float64 = 0.0, ht_ml::Float64 = 0.0, ht_mh::Float64 = 0.0, ht_ef::Float64 = 0.0, 
        ht_el::Float64 = 0.0, ht_eh::Float64 = 0.0, 
        
        rhot_ps::Float64 = 0.0, rhot_pt::Float64 = 0.0, rhot_ss::Float64 = 0.0, rhot_st::Float64 = 0.0, 
        rhot_Ts::Float64 = 0.0, rhot_Tt::Float64 = 0.0, rhot_hs::Float64 = 0.0, rhot_ht::Float64 = 0.0, 
        rhot_Mi::Float64 = 0.0, rhot_pf::Float64 = 0.0, rhot_pl::Float64 = 0.0, rhot_ph::Float64 = 0.0, 
        rhot_mf::Float64 = 0.0, rhot_ml::Float64 = 0.0, rhot_mh::Float64 = 0.0, rhot_ef::Float64 = 0.0, 
        rhot_el::Float64 = 0.0, rhot_eh::Float64 = 0.0, 
        
        st_ps::Float64 = 0.0, st_pt::Float64 = 0.0, st_ss::Float64 = 0.0, st_st::Float64 = 0.0, 
        st_Ts::Float64 = 0.0, st_Tt::Float64 = 0.0, st_hs::Float64 = 0.0, st_ht::Float64 = 0.0, 
        st_Mi::Float64 = 0.0, st_pf::Float64 = 0.0, st_pl::Float64 = 0.0, st_ph::Float64 = 0.0, 
        st_mf::Float64 = 0.0, st_ml::Float64 = 0.0, st_mh::Float64 = 0.0, st_ef::Float64 = 0.0, 
        st_el::Float64 = 0.0, st_eh::Float64 = 0.0, 

    )
        return new(
            # ---- Physical Flow
            Mi, A, u, u_Mi,

            # ---- Static Quantities
            ps, Ts, hs, rhos, ss, cps, cvs, gams, as,

            # ---- Total Quantities
            pt, Tt, ht, rhot, st, cpt, cvt, gamt, at,
            
            # ---- Static Derivatives
            ps_ps, ps_pt, ps_ss, ps_st, ps_Ts, ps_Tt, ps_hs, ps_ht, ps_Mi, ps_pf, ps_pl, ps_ph, ps_mf, ps_ml, ps_mh, ps_ef, ps_el, ps_eh, 
            Ts_ps, Ts_pt, Ts_ss, Ts_st, Ts_Ts, Ts_Tt, Ts_hs, Ts_ht, Ts_Mi, Ts_pf, Ts_pl, Ts_ph, Ts_mf, Ts_ml, Ts_mh, Ts_ef, Ts_el, Ts_eh, 
            hs_ps, hs_pt, hs_ss, hs_st, hs_Ts, hs_Tt, hs_hs, hs_ht, hs_Mi, hs_pf, hs_pl, hs_ph, hs_mf, hs_ml, hs_mh, hs_ef, hs_el, hs_eh, 
            rhos_ps, rhos_pt, rhos_ss, rhos_st, rhos_Ts, rhos_Tt, rhos_hs, rhos_ht, rhos_Mi, rhos_pf, rhos_pl, rhos_ph, rhos_mf, rhos_ml, rhos_mh, rhos_ef, rhos_el, rhos_eh, 
            ss_ps, ss_pt, ss_ss, ss_st, ss_Ts, ss_Tt, ss_hs, ss_ht, ss_Mi, ss_pf, ss_pl, ss_ph, ss_mf, ss_ml, ss_mh, ss_ef, ss_el, ss_eh, 
            
            # ---- Total Derivatives
            pt_ps, pt_pt, pt_ss, pt_st, pt_Ts, pt_Tt, pt_hs, pt_ht, pt_Mi, pt_pf, pt_pl, pt_ph, pt_mf, pt_ml, pt_mh, pt_ef, pt_el, pt_eh, 
            Tt_ps, Tt_pt, Tt_ss, Tt_st, Tt_Ts, Tt_Tt, Tt_hs, Tt_ht, Tt_Mi, Tt_pf, Tt_pl, Tt_ph, Tt_mf, Tt_ml, Tt_mh, Tt_ef, Tt_el, Tt_eh, 
            ht_ps, ht_pt, ht_ss, ht_st, ht_Ts, ht_Tt, ht_hs, ht_ht, ht_Mi, ht_pf, ht_pl, ht_ph, ht_mf, ht_ml, ht_mh, ht_ef, ht_el, ht_eh, 
            rhot_ps, rhot_pt, rhot_ss, rhot_st, rhot_Ts, rhot_Tt, rhot_hs, rhot_ht, rhot_Mi, rhot_pf, rhot_pl, rhot_ph, rhot_mf, rhot_ml, rhot_mh, rhot_ef, rhot_el, rhot_eh, 
            st_ps, st_pt, st_ss, st_st, st_Ts, st_Tt, st_hs, st_ht, st_Mi, st_pf, st_pl, st_ph, st_mf, st_ml, st_mh, st_ef, st_el, st_eh, 
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
end # turbine