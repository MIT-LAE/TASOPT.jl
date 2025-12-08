@testset "mission deviation" verbose=true begin
    default_ac = load_default_model()
    size_aircraft!(default_ac; printiter = false)

    mission_idx = 1
    cruise_nodes = (
        TASOPT.ipcruise1,
        TASOPT.ipdeviation1,
        TASOPT.ipdeviation2,
        TASOPT.ipdeviation3,
        TASOPT.ipdeviation4,
        TASOPT.ipcruise2,
    )

    function get_deviation_info(ac, imission)
        range_start = ac.para[TASOPT.iaRange, TASOPT.ipcruise1, imission]
        alt_start = ac.para[TASOPT.iaalt, TASOPT.ipcruise1, imission]
        gam_base = ac.para[TASOPT.iagamV, TASOPT.ipcruise1, imission]

        s_vals = [
            ac.para[TASOPT.iaRange, ip, imission] - range_start for ip in cruise_nodes
        ]
        offsets = [
            ac.para[TASOPT.iaalt, ip, imission] - (alt_start + gam_base * s)
            for (ip, s) in zip(cruise_nodes, s_vals)
        ]

        gam_up = ac.para[TASOPT.iagamV, TASOPT.ipdeviation1, imission]
        gam_dn = ac.para[TASOPT.iagamV, TASOPT.ipdeviation3, imission]

        base_length = (s_vals[3] - s_vals[2]) + (s_vals[5] - s_vals[4])
        R_hold = s_vals[4] - s_vals[3]

        return (
            s_vals = s_vals,
            offsets = offsets,
            gam_base = gam_base,
            gam_up = gam_up,
            gam_dn = gam_dn,
            base_length = base_length,
            R_hold = R_hold,
            dRcruise = s_vals[end],
            Δh_eff = offsets[3],
        )
    end

    function gamma_cap(ac, imission)
        WMTO = ac.parg[TASOPT.igWMTO]
        BW_cap =
            ac.para[TASOPT.iafracW, TASOPT.ipcruise1, imission] * WMTO +
            ac.para[TASOPT.iaWbuoy, TASOPT.ipcruise1, imission]
        DoL_cap =
            ac.para[TASOPT.iaCD, TASOPT.ipcruise1, imission] /
            ac.para[TASOPT.iaCL, TASOPT.ipcruise1, imission]
        F_available =
            ac.pare[TASOPT.ieFe, TASOPT.ipclimbn, imission] * ac.parg[TASOPT.igneng]
        sin_cap = F_available / BW_cap - DoL_cap
        return asin(clamp(sin_cap, -1.0, 1.0))
    end

    @testset "Sample deviation test" begin
        ac_default = deepcopy(default_ac)

        ac_default.parm[TASOPT.imDeviationHeight, mission_idx] = 304.8 # 1000 ft
        ac_default.parm[TASOPT.imDeviationLength, mission_idx] = 185200 # 100 nmi
        ac_default.parm[TASOPT.imDeviationStartFromTOC, mission_idx] = 185200 # 100 nmi

        TASOPT._mission_iteration!(ac_default, mission_idx, false)

        snap = get_deviation_info(ac_default, mission_idx)

        Δh_req = ac_default.parm[TASOPT.imDeviationHeight, mission_idx]
        R_dev_req = ac_default.parm[TASOPT.imDeviationLength, mission_idx]
        S_dev_req = ac_default.parm[TASOPT.imDeviationStartFromTOC, mission_idx]

        gam_base = snap.gam_base
        gam_up = snap.gam_up
        gam_dn = snap.gam_dn
        den_up = gam_up - gam_base
        den_dn = gam_dn - gam_base

        Δs_up = snap.s_vals[3] - snap.s_vals[2]
        Δs_dn = snap.s_vals[5] - snap.s_vals[4]
        base_length = snap.base_length
        R_hold = snap.R_hold
        dR_available = snap.dRcruise
        R_dev_total = base_length + R_hold
        S_dev_result = snap.s_vals[2]
        S_dev_max = max(dR_available - R_dev_total, 0.0)

        denom = abs((1 / den_up) - (1 / den_dn))
        Δh_limit = dR_available / denom
        Δh_eff = clamp(Δh_req, -Δh_limit, Δh_limit)

        γ_climb = ac_default.para[TASOPT.iagamV, TASOPT.ipclimbn, mission_idx]
        γ_eff = min(abs(γ_climb), gamma_cap(ac_default, mission_idx))
        sign_h = Δh_eff == 0.0 ? 0.0 : sign(Δh_eff)
        expected_up = sign_h * γ_eff
        expected_dn = -sign_h * γ_eff

        @test isapprox(snap.offsets[1], 0.0; atol = 1e-6)
        @test isapprox(snap.offsets[2], 0.0; atol = 1e-6)
        @test isapprox(snap.offsets[3], Δh_eff; atol = 1e-6, rtol = 1e-6)
        @test isapprox(snap.offsets[4], Δh_eff; atol = 1e-6, rtol = 1e-6)
        @test isapprox(snap.offsets[5], 0.0; atol = 1e-6)
        @test isapprox(snap.offsets[6], 0.0; atol = 1e-6)

        @test isapprox(S_dev_result, min(S_dev_req, S_dev_max); atol = 1e-6)
        @test isapprox(R_hold, max(R_dev_req - base_length, 0.0); atol = 1e-6)

        @test isapprox(Δs_up * den_up, Δh_eff; atol = 1e-6, rtol = 1e-6)
        @test isapprox(Δs_dn * den_dn, -Δh_eff; atol = 1e-6, rtol = 1e-6)
        @test isapprox(snap.Δh_eff, Δh_eff; atol = 1e-6)

        @test isapprox(gam_up, expected_up; atol = 1e-8)
        @test isapprox(gam_dn, expected_dn; atol = 1e-8)
    end

    @testset "No deviation (default model)" begin
        ac_default = deepcopy(default_ac)

        TASOPT._mission_iteration!(ac_default, mission_idx, false)

        snap = get_deviation_info(ac_default, mission_idx)

        Δh_req = ac_default.parm[TASOPT.imDeviationHeight, mission_idx]
        R_dev_req = ac_default.parm[TASOPT.imDeviationLength, mission_idx]
        S_dev_req = ac_default.parm[TASOPT.imDeviationStartFromTOC, mission_idx]

        gam_base = snap.gam_base
        gam_up = snap.gam_up
        gam_dn = snap.gam_dn
        den_up = gam_up - gam_base
        den_dn = gam_dn - gam_base

        Δs_up = snap.s_vals[3] - snap.s_vals[2]
        Δs_dn = snap.s_vals[5] - snap.s_vals[4]
        base_length = snap.base_length
        R_hold = snap.R_hold
        dR_available = snap.dRcruise
        R_dev_total = base_length + R_hold
        S_dev_result = snap.s_vals[2]
        S_dev_max = max(dR_available - R_dev_total, 0.0)

        denom = abs((1 / den_up) - (1 / den_dn))
        Δh_limit = dR_available / denom
        Δh_eff = clamp(Δh_req, -Δh_limit, Δh_limit)

        γ_cruise = ac_default.para[TASOPT.iagamV, TASOPT.ipcruise1, mission_idx]
        γ_eff = min(abs(γ_cruise), gamma_cap(ac_default, mission_idx))
        if Δh_eff == 0.0
            expected_up = γ_eff
            expected_dn = γ_eff
        else
            sign_h = sign(Δh_eff)
            expected_up = sign_h * γ_eff
            expected_dn = -sign_h * γ_eff
        end
        @test isapprox(snap.offsets[1], 0.0; atol = 1e-6)
        @test isapprox(snap.offsets[2], 0.0; atol = 1e-6)
        @test isapprox(snap.offsets[3], Δh_eff; atol = 1e-6, rtol = 1e-6)
        @test isapprox(snap.offsets[4], Δh_eff; atol = 1e-6, rtol = 1e-6)
        @test isapprox(snap.offsets[5], 0.0; atol = 1e-6)
        @test isapprox(snap.offsets[6], 0.0; atol = 1e-6)

        @test isapprox(S_dev_result, R_hold; atol = 1e-6)

        @test isapprox(Δs_up * den_up, Δh_eff; atol = 1e-6, rtol = 1e-6)
        @test isapprox(Δs_dn * den_dn, -Δh_eff; atol = 1e-6, rtol = 1e-6)
        @test isapprox(snap.Δh_eff, Δh_eff; atol = 1e-6)

        @test isapprox(gam_up, expected_up; atol = 1e-8)
        @test isapprox(gam_dn, expected_dn; atol = 1e-8)
    end

    @testset "Deviation start clamps to cruise span" begin
        ac_start = deepcopy(default_ac)
        ac_start.parm[TASOPT.imDeviationHeight, mission_idx] = 304.8 # ensure deviation active
        ac_start.parm[TASOPT.imDeviationLength, mission_idx] = 185200.0
        ac_start.parm[TASOPT.imDeviationStartFromTOC, mission_idx] = 1.0e7

        TASOPT._mission_iteration!(ac_start, mission_idx, false)
        snap = get_deviation_info(ac_start, mission_idx)

        S_dev_req = ac_start.parm[TASOPT.imDeviationStartFromTOC, mission_idx]
        R_dev_req = ac_start.parm[TASOPT.imDeviationLength, mission_idx]
        Δh_req = ac_start.parm[TASOPT.imDeviationHeight, mission_idx]

        gam_base = snap.gam_base
        gam_up = snap.gam_up
        gam_dn = snap.gam_dn
        den_up = gam_up - gam_base
        den_dn = gam_dn - gam_base

        base_length = snap.base_length
        R_hold = snap.R_hold
        dR_available = snap.dRcruise
        R_dev_total = base_length + R_hold
        S_dev_max = max(dR_available - R_dev_total, 0.0)

        @test S_dev_req > S_dev_max
        @test isapprox(snap.s_vals[2], S_dev_max; atol = 1e-6)

        denom = abs((1 / den_up) - (1 / den_dn))
        Δh_limit = dR_available / denom
        Δh_eff = clamp(Δh_req, -Δh_limit, Δh_limit)

        @test isapprox(snap.offsets[1], 0.0; atol = 1e-6)
        @test isapprox(snap.offsets[6], 0.0; atol = 1e-6)
        @test isapprox(snap.offsets[3], Δh_eff; atol = 1e-6, rtol = 1e-6)
        @test isapprox(snap.offsets[4], Δh_eff; atol = 1e-6, rtol = 1e-6)

        @test isapprox(R_hold, max(R_dev_req - base_length, 0.0); atol = 1e-6)
    end

    @testset "Altitude change clamps to cruise span limit" begin
        ac_clamp = deepcopy(default_ac)
        ac_clamp.parm[TASOPT.imDeviationLength, mission_idx] = 185200 # 100 nmi
        ac_clamp.parm[TASOPT.imDeviationStartFromTOC, mission_idx] = 185200 # 100 nmi
        requested = 3.0e4
        ac_clamp.parm[TASOPT.imDeviationHeight, mission_idx] = requested

        TASOPT._mission_iteration!(ac_clamp, mission_idx, false)
        snap = get_deviation_info(ac_clamp, mission_idx)

        gam_base = snap.gam_base
        gam_up = snap.gam_up
        gam_dn = snap.gam_dn
        den_up = gam_up - gam_base
        den_dn = gam_dn - gam_base

        denom = abs((1 / den_up) - (1 / den_dn))
        Δh_limit = snap.dRcruise / denom
        abs_offset = abs(snap.offsets[3])
        expected_abs = min(abs(requested), Δh_limit)
        @test isapprox(abs_offset, expected_abs; atol = 1e-6, rtol = 1e-6)
        @test snap.offsets[3] > 0.0
        @test isapprox(snap.offsets[3], snap.offsets[4]; atol = 1e-6)

        γ_climb = ac_clamp.para[TASOPT.iagamV, TASOPT.ipclimbn, mission_idx]
        γ_eff = min(abs(γ_climb), gamma_cap(ac_clamp, mission_idx))
        expected_up = γ_eff

        @test isapprox(gam_up, expected_up; atol = 1e-8)
        @test isapprox(gam_dn, -expected_up; atol = 1e-8)
    end

    @testset "Hold segment disappears when requested length is too short" begin
        ac_short = deepcopy(default_ac)
        ac_short.parm[TASOPT.imDeviationHeight, mission_idx] = 304.8
        ac_short.parm[TASOPT.imDeviationLength, mission_idx] = 1852.0 * 1
        ac_short.parm[TASOPT.imDeviationStartFromTOC, mission_idx] = 1.0e4

        TASOPT._mission_iteration!(ac_short, mission_idx, true)
        snap = get_deviation_info(ac_short, mission_idx)

        @test snap.base_length > 0.0
        @test isapprox(snap.R_hold, 0.0; atol = 1e-6)
        @test isapprox(snap.s_vals[4], snap.s_vals[3]; atol = 1e-6)

        S_dev_req = ac_short.parm[TASOPT.imDeviationStartFromTOC, mission_idx]
        gam_base = snap.gam_base
        gam_up = snap.gam_up
        gam_dn = snap.gam_dn
        den_up = gam_up - gam_base
        den_dn = gam_dn - gam_base

        denom = abs((1 / den_up) - (1 / den_dn))
        Δh_limit = snap.dRcruise / denom
        Δh_req = ac_short.parm[TASOPT.imDeviationHeight, mission_idx]
        Δh_eff = clamp(Δh_req, -Δh_limit, Δh_limit)

        S_dev_max = max(snap.dRcruise - snap.base_length, 0.0)
        @test isapprox(snap.s_vals[2], clamp(S_dev_req, 0.0, S_dev_max); atol = 1e-6)
        @test isapprox(snap.offsets[3], Δh_eff; atol = 1e-6, rtol = 1e-6)
    end
end
