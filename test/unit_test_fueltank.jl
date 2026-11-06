#Sample parameters
ac = load_default_model()

ifuel = 40
z = 11e3
Mair = 0.8
xftank = 15.0
ac.options.ifuel = ifuel
ac.para[iaalt, ipcruise1,1] = z
ac.para[iaMach, ipcruise1,1] = Mair
ac.parg[igxftank] = xftank
ac.fuse_tank.tank_count = 1
ac.parg[igWfuel] = 1e5

fuse_tank = ac.fuse_tank
fuse_tank.fueltype = "LH2"
fuse_tank.placement = "front"

fuse_tank.clearance_fuse = 0.1
TSL = 288.2
fuse_tank.TSLtank = [TSL]

fuse_tank.pvent = 2e5
fuse_tank.ARtank = 2.0
fuse_tank.inner_material = TASOPT.StructuralAlloy("Al-2219-T87")

fuse_tank.t_insul = [0.15,0.15,0.15]
fuse_tank.material_insul = [ThermalInsulator("rohacell41s"), ThermalInsulator("rohacell41s"), ThermalInsulator("polyurethane27")]
fuse_tank.sizes_insulation = true
fuse_tank.iinsuldes = [1,2,3]
fuse_tank.boiloff_rate = 0.15

fuse_tank.ullage_frac = 0.1
fuse_tank.qfac = 1.3
fuse_tank.theta_inner = 1.0
fuse_tank.ftankadd = 0.1
fuse_tank.ew = 0.9

fuse_tank.outer_material = StructuralAlloy("Al-2219-T87")
fuse_tank.theta_outer = [1.0, 2.0]
fuse_tank.Ninterm = 1.0

#Calculate fuel temperature and density as a function of pressure
β0 = 1 - fuse_tank.ullage_frac
fuel_mix = TASOPT.SaturatedMixture(fuse_tank.fueltype, fuse_tank.pvent, β0)

Tsat = fuel_mix.liquid.Tsat
ρl = fuel_mix.liquid.ρ
ρg = fuel_mix.gas.ρ
hvap = fuel_mix.hvap

fuse_tank.rhofuel = ρl
fuse_tank.rhofuelgas = ρg
fuse_tank.Tfuel = Tsat
fuse_tank.hvap = hvap

#Fuselage
fuse = ac.fuselage
fuse.layout.cross_section.radius = 2.5
fuse.layout.cross_section.bubble_lower_downward_shift = 0.3

@testset "Fuselage tank" begin
    @testset "Foam insulation" begin

        TASOPT.tanksize!(ac, 1)
        inner_tank_mech = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)

        @test inner_tank_mech.Wtank       ≈ 60109.124334566346
        @test inner_tank_mech.Winsul_sum  ≈ 38147.81559388217
        @test inner_tank_mech.Vfuel       ≈ 166.77327116515787
        @test inner_tank_mech.Shead_insul ≈ [16.42133035129676, 20.311090752967278, 24.620380536430755, 29.34797626786876]
        @test inner_tank_mech.Rtank_outer ≈ 1.8743852420176998
        @test inner_tank_mech.l_tank      ≈ 15.771807531012154
        @test inner_tank_mech.l_cyl       ≈ 12.842701653674803
        
        outputs_thermal = TASOPT.CryoTank.tankWthermal(fuse, fuse_tank, z, TSL, Mair, xftank, ifuel)

        outputs_thermal_check = (1835.8608916011449,)

        for i in 1:length(outputs_thermal)
            @test outputs_thermal[i] ≈ outputs_thermal_check[i]
        end
    end
    
    fuse_tank.t_insul = [0.065]
    fuse_tank.material_insul = [ThermalInsulator("vacuum")]
    fuse_tank.sizes_insulation = false

    @testset "Vacuum insulation" begin
        TASOPT.tanksize!(ac, 1)

        inner_tank_vac = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)

        @test inner_tank_vac.Wtank       ≈ 24229.845334360776
        @test inner_tank_vac.Winsul_sum  ≈ 0.0
        @test inner_tank_vac.Vfuel       ≈ 166.77327116515787
        @test inner_tank_vac.Shead_insul ≈ [25.483812169139977, 27.21328988831087]
        @test inner_tank_vac.Rtank_outer ≈ 2.335
        @test inner_tank_vac.l_tank      ≈ 9.994915110929762
        @test inner_tank_vac.l_cyl       ≈ 7.525566077704729

        Winnertank = inner_tank_vac.Wtank
        l_cyl      = inner_tank_vac.l_cyl
        Winner_tot = Winnertank + fuse_tank.Wfuelintank

        fuse_tank.Ninterm = 1.0
        Ninterm = TASOPT.CryoTank.optimize_outer_tank(fuse, fuse_tank, Winner_tot, l_cyl)
        Ninterm_check = 14.025390625
        @test Ninterm ≈ Ninterm_check

        outer_tank_vac = TASOPT.CryoTank.size_outer_tank(fuse, fuse_tank, Winner_tot, l_cyl, Ninterm_check)
        @test outer_tank_vac.Wtank   ≈ 87481.40843690369
        @test outer_tank_vac.Wcyl    ≈ 31381.827647692957
        @test outer_tank_vac.Whead   ≈ 13786.613724211644
        @test outer_tank_vac.Wstiff  ≈ 20573.498028341648
        @test outer_tank_vac.Souter  ≈ 171.6223477616751
        @test outer_tank_vac.Shead   ≈ 26.90237940128641
        @test outer_tank_vac.Scyl    ≈ 117.81758895910225
        @test outer_tank_vac.t_cyl   ≈ 0.009560503271690743
        @test outer_tank_vac.t_head  ≈ 0.018394143361195915
        @test outer_tank_vac.l_outer ≈ 9.96235436442712

    end
    
    @testset "Stringer sizing" begin
        θ = 85 * pi /180
        outputs_bendM = TASOPT.CryoTank.stiffeners_bendingM(θ)
        outputs_bendM_check = [1.1749999999999998, 0.09494892091249901]
        for i in 1:length(outputs_bendM)
            @test outputs_bendM[i] ≈ outputs_bendM_check[i]
        end

        θ1 = 45 * pi /180
        θ2 = 135 * pi /180
        outputs_bendM_outer = TASOPT.CryoTank.stiffeners_bendingM_outer(θ1, θ2)
        outputs_bendM_outer_check = [pi/2, 0.6166583566771067]
        for i in 1:length(outputs_bendM)
            @test outputs_bendM_outer[i] ≈ outputs_bendM_outer_check[i]
        end

        perim_vessel = fuse.layout.cross_section.perimeter

        Wstiff = TASOPT.CryoTank.stiffener_weight("outer", 7e4, fuse.layout.radius, perim_vessel, fuse_tank.outer_material.UTS / 4, 
        fuse_tank.outer_material.ρ, θ1, θ2, 10.0, 5.0, fuse_tank.outer_material.E)
        Wstiff_check = 1717.0752091513864
        
        @test Wstiff ≈ Wstiff_check
    end

    @testset "Thermal models" begin
        outputs_h = TASOPT.CryoTank.tank_heat_coeff(21.0, ifuel, fuse_tank.Tfuel, 5.0)
        outputs_h_check = (133.43555789562762)
        for i in 1:length(outputs_h)
            @test outputs_h[i] ≈ outputs_h_check[i]
        end

        outputs_free = TASOPT.CryoTank.freestream_heat_coeff(z, TSL, Mair, xftank, 240.0)
        outputs_free_check = (91.30101055991513, 218.06145060705106, 243.26644796151263)
        for i in 1:length(outputs_free)
            @test outputs_free[i] ≈ outputs_free_check[i]
        end

        outputs_natural = TASOPT.CryoTank.freestream_heat_coeff(0.0, TSL, 0.0, xftank, 287.0, fuse.layout.radius)
        outputs_natural_check = (1.6795624846472137, 288.2, 288.2)
        for i in 1:length(outputs_natural)
            @test outputs_natural[i] ≈ outputs_natural_check[i]
        end

        Taw = outputs_free_check[3]
        Rvac = TASOPT.CryoTank.vacuum_resistance(fuse_tank.Tfuel, Taw, 90.0, 100.0)
        Rvac_check = 0.3565313051181731

        @test Rvac ≈ Rvac_check
    end
end

@testset "Tank sizing lower bounds" begin
    CT = TASOPT.CryoTank
    mat = TASOPT.materials.StructuralAlloy("Al-2219-T87")
    s_a = mat.UTS / 4
    θ1, θ2 = 45 * pi / 180, 135 * pi / 180
    perim = 2π * 1.86

    @testset "stiffener_weight light load" begin
        Ws = [CT.stiffener_weight("inner", W, 1.86, perim, s_a, mat.ρ, 1.0) for W in (0.0, 1e3, 5e3)]
        @test all(isfinite, Ws)
        @test all(>(0), Ws)
        @test Ws[1] ≈ Ws[2] ≈ Ws[3] #Constant minimum beam below the threshold
        @test CT.stiffener_weight("inner", 1e6, 1.86, perim, s_a, mat.ρ, 1.0) > Ws[1] #Monotonic in load
        @test CT.stiffener_weight("outer", 0.0, 1.86, perim, s_a, mat.ρ, θ1, θ2, 3.0, 1.0, mat.E) > 0
    end

    @testset "outer_cylinder_thickness_ratio" begin
        pc = 4 * TASOPT.pref
        ν = mat.ν
        E = mat.E
        prev = 0.0
        for L_Do in (1e-3, 1e-2, 0.1, 0.5, 1.0)
            t_Do = CT.outer_cylinder_thickness_ratio(E, ν, L_Do, pc)
            @test 0 < t_Do < (L_Do / 0.45)^2
            res = 2.42 * E * t_Do^(5/2) / ((1 - ν^2)^(3/4) * (L_Do - 0.45 * sqrt(t_Do))) - pc
            @test abs(res) < 1e-3 * pc
            @test t_Do > prev #Monotonic in L/D
            prev = t_Do
        end
        @test_throws DomainError CT.outer_cylinder_thickness_ratio(E, ν, 0.0, pc)
        @test_throws DomainError CT.outer_cylinder_thickness_ratio(E, ν, -1.0, pc)
    end

    @testset "insulation thickness checks" begin
        @test CT.check_insulation_thickness([0.0, 0.1]) === nothing
        @test_throws DomainError CT.check_insulation_thickness([0.1, -1e-3])
        @test_throws DomainError CT.check_insulation_thickness([0.1, NaN])
        t = [-0.1, 0.05, 0.0]
        @test_logs (:warn,) CT.clamp_insulation_thickness!(t)
        @test t == [CT.T_INSUL_MIN, 0.05, CT.T_INSUL_MIN]
    end

    @testset "residuals_Q: non-positive temperature with a vacuum layer" begin
        #A Newton step can propose a non-positive interface temperature, for which the vacuum resistance has no real value. 
        #The residual must then return the large-residual sentinel instead of throwing (throwing is ~100x slower).
        #Uses the shared fixture; the fields that are changed are restored afterwards
        ft, fs = fuse_tank, fuse
        saved = (ft.material_insul, ft.t_insul, ft.sizes_insulation, ft.iinsuldes)
        ft.material_insul = [ThermalInsulator("vacuum"), ThermalInsulator("polyurethane32")]
        ft.t_insul = [0.02, 0.05]
        ft.iinsuldes = [1, 2]
        ft.sizes_insulation = false
        inner = CT.size_inner_tank(fs, ft, ft.t_insul)
        p = CT.thermal_params{typeof(fs.layout.cross_section), Float64}()
        p.l_cyl = inner.l_cyl; p.l_tank = inner.l_tank; p.r_tank = inner.Rtank_outer; p.Shead = inner.Shead_insul
        p.t_cond = ft.t_insul; p.material = ft.material_insul; p.Tfuel = ft.Tfuel
        p.z = 11e3; p.TSL = TSL; p.Mair = 0.8; p.xftank = 15.0; p.ifuel = ifuel; p.fuse_cs = fs.layout.cross_section

        x_bad = [5000.0, 21.0, -10.0, 250.0] #Q, T_wall, T_ins (vacuum outer face < 0), T_fuselage
        @test (@test_nowarn CT.residuals_Q(x_bad, p, "Q_unknown")) == fill(1e3, 4)
        x_zero = [5000.0, 21.0, 0.0, 250.0] #Boundary: exactly 0 K is also not a valid vacuum-side temperature
        @test CT.residuals_Q(x_zero, p, "Q_unknown") == fill(1e3, 4)
        p.Q = 1000.0
        @test CT.residuals_Q([21.0, -10.0, 250.0], p, "Q_known") == fill(1e3, 4) #Q_known mode has one more residual than states
        x_ok = [500.0, 21.0, 100.0, 250.0]
        r = CT.residuals_Q(x_ok, p, "Q_unknown")
        @test all(isfinite, r) && !all(==(1e3), r)

        ft.material_insul, ft.t_insul, ft.sizes_insulation, ft.iinsuldes = saved
    end
end

@testset "Continuous tank sizing" begin
    CT = TASOPT.CryoTank
    fuse_tank.t_insul = [0.065]
    fuse_tank.material_insul = [ThermalInsulator("vacuum")]
    fuse_tank.sizes_insulation = false
    Rmax = fuse.layout.radius - sum(fuse_tank.t_insul) - fuse_tank.clearance_fuse
    Wf_orig = fuse_tank.Wfuelintank

    #Sweep the fuel weight over the sphere, morphing and cylinder regimes
    Wfs = 10 .^ range(1, 6, length = 400) .* gee
    geoms = map(Wfs) do Wf
        fuse_tank.Wfuelintank = Wf
        CT.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
    end
    fuse_tank.Wfuelintank = Wf_orig

    R  = [g.Rtank_outer for g in geoms]
    L  = [g.l_cyl for g in geoms]
    AR = [g.AR_head for g in geoms]
    Wt = [g.Wtank for g in geoms]
    lt = [g.l_tank for g in geoms]

    @test all(L .>= 0)
    @test all(1 .<= AR .<= fuse_tank.ARtank + 1e-12)
    @test all(R .<= Rmax + 1e-12)
    @test R[1] < Rmax && R[end] ≈ Rmax #Starts as a small sphere and ends at the maximum radius
    @test L[1] == 0.0 && AR[1] == 1.0 #Spherical at low volume
    @test AR[end] ≈ fuse_tank.ARtank && L[end] > 0 #Design heads plus cylinder at high volume

    #Monotonic growth and no jumps (largest relative step between consecutive points is small)
    for v in (R, Wt, lt)
        @test issorted(v)
        @test maximum(diff(v) ./ v[1:end-1]) < 0.1 #Fuel weight steps are about 3 % apart
    end
    @test issorted(L) && issorted(AR)
    @test maximum(abs.(diff(AR))) < 0.1
    @test maximum(diff(L)[L[1:end-1] .> 0]) < 0.2 * maximum(L) 

    #Weights reproduce the original (cylinder-regime) results where the old model applied
    fuse_tank.Wfuelintank = Wf_orig
    g = CT.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
    @test g.AR_head == fuse_tank.ARtank
    @test g.Rtank_outer ≈ Rmax

    #Full tank sizing with an outer vessel stays finite and continuous across regimes
    fuse_tank.sizes_insulation = false
    Wfuels = 10 .^ range(3, 5.6, length = 300) #about 2 % apart
    Wtanks = Float64[]
    for Wfuel in Wfuels
        ac.parg[igWfuel] = Wfuel
        TASOPT.tanksize!(ac, 1)
        @test isfinite(ac.parg[igWftank]) && ac.parg[igWftank] > 0
        @test isfinite(ac.parg[iglftank]) && ac.parg[iglftank] > 0
        push!(Wtanks, ac.parg[igWftank])
    end
    @test issorted(Wtanks)
    @test maximum(diff(Wtanks) ./ Wtanks[1:end-1]) < 0.1 #No jumps in tank weight between regimes
    ac.parg[igWfuel] = 1e5
end


@testset "Unreachable boil-off target" begin
    CT = TASOPT.CryoTank
    #Single vacuum layer: its heat leak is limited, so a very high boil-off target is unreachable
    fuse_tank.material_insul = [ThermalInsulator("vacuum")]
    fuse_tank.t_insul = [0.065]
    fuse_tank.iinsuldes = [1]
    fuse_tank.sizes_insulation = true
    fuse_tank.last_MLI_sol = Float64[]
    fuse_tank.last_thermal_sol = Float64[]
    ac.parg[igWfuel] = 1e5

    zc = ac.para[iaalt, ipcruise1, 1]
    Mc = ac.para[iaMach, ipcruise1, 1]

    #Highest attainable heat rate is consistent with the full thermal model for the same stack
    fuse_tank.boiloff_rate = 1e3
    fuse_tank.sizes_insulation = false
    TASOPT.tanksize!(ac, 1) #Stores the tank geometry and position used by tankWthermal
    xc = ac.parg[igxftank]
    Q_max = CT.max_heat_rate(fuse, fuse_tank, 0.0, zc, TSL, Mc, xc, ifuel)
    @test isfinite(Q_max) && Q_max > 0
    @test CT.tankWthermal(fuse, fuse_tank, zc, TSL, Mc, xc, ifuel) ≈ fuse_tank.qfac * Q_max rtol = 1e-3

    #Unreachable target: no error, warns, and stack is physical
    fuse_tank.t_insul = [0.065]
    fuse_tank.sizes_insulation = true
    fuse_tank.boiloff_rate = 1e3 #%/h
    @test CT.boiloff_heat_rate(fuse_tank) > fuse_tank.qfac * Q_max
    @test_logs (:warn, r"cannot be reached") match_mode = :any TASOPT.tanksize!(ac, 1)
    @test all(fuse_tank.t_insul .>= CT.T_INSUL_MIN)
    @test fuse_tank.t_insul[1] ≈ CT.T_INSUL_MIN
    @test isfinite(ac.parg[igWftank]) && ac.parg[igWftank] > 0

    #Reachable target: normal solve (no warning about reachability), thickness above the minimum
    fuse_tank.t_insul = [0.065]
    fuse_tank.material_insul = [ThermalInsulator("rohacell41s"), ThermalInsulator("rohacell41s"), ThermalInsulator("polyurethane27")]
    fuse_tank.t_insul = [0.15, 0.15, 0.15]
    fuse_tank.iinsuldes = [1, 2, 3]
    fuse_tank.boiloff_rate = 0.15
    fuse_tank.last_MLI_sol = Float64[]
    fuse_tank.last_thermal_sol = Float64[]
    TASOPT.tanksize!(ac, 1)
    @test all(fuse_tank.t_insul .> 10 * CT.T_INSUL_MIN)
end

@testset "Thermal solve robustness" begin
    CT = TASOPT.CryoTank
    ft = fuse_tank
    zc = ac.para[iaalt, ipcruise1, 1]
    Mc = ac.para[iaMach, ipcruise1, 1]
    ac.parg[igWfuel] = 1e5

    function setup_stack!(materials, t)
        ft.material_insul = ThermalInsulator.(materials)
        ft.t_insul = copy(t)
        ft.iinsuldes = collect(1:length(t))
        ft.sizes_insulation = false
        ft.last_MLI_sol = Float64[]
        ft.last_thermal_sol = Float64[]
        TASOPT.tanksize!(ac, 1) #Stores the tank geometry and position used by the thermal model
        return ac.parg[igxftank]
    end

    @testset "max_heat_rate: thinnest admissible non-vacuum stack" begin
        #The feasibility pre-check evaluates the stack reduced to T_INSUL_MIN, whose resistance is orders of magnitude below that of the
        #design stack. The solve must still find the physical root (the pre-check returning NaN silently disables it).
        xc = setup_stack!(["polyurethane32", "polyurethane32", "polyurethane32"], [0.15, 0.15, 0.15])
        Δt_min = CT.min_insulation_shift(ft.t_insul, ft.iinsuldes)
        @test Δt_min ≈ CT.T_INSUL_MIN - 0.15
        Q_thin = CT.max_heat_rate(fuse, ft, Δt_min, zc, TSL, Mc, xc, ifuel)
        Q_full = CT.max_heat_rate(fuse, ft, 0.0, zc, TSL, Mc, xc, ifuel)
        @test isfinite(Q_thin) && isfinite(Q_full)
        @test Q_thin > 10 * Q_full > 0 #Much thinner insulation leaks much more
    end

    @testset "tankWthermal: non-physical warm start is not accepted" begin
        #With a vacuum layer the thermal network has non-physical stagnation points (Q ~ 0 and a wall below the fuel temperature).
        #A stale warm start there must not be returned as a solution, or stored for the next call.
        xc = setup_stack!(["vacuum", "polyurethane32", "polyurethane32"], [0.02, 0.1, 0.1])
        Q_ref = CT.tankWthermal(fuse, ft, zc, TSL, Mc, xc, ifuel)
        @test isfinite(Q_ref) && Q_ref > 0

        _, _, Taw = CT.freestream_heat_coeff(zc, TSL, Mc, xc)
        ft.last_thermal_sol = [3731.0, 11.9, 1e-19, 152.0, 240.0] #Non-physical stagnation point (wall below the fuel temperature), as seen in practice
        Q_poisoned = CT.tankWthermal(fuse, ft, zc, TSL, Mc, xc, ifuel)
        @test Q_poisoned ≈ Q_ref rtol = 1e-3
        @test ft.last_thermal_sol[1] > 1.0 #Cached state is a true root, not the stagnation point
        @test all(T -> ft.Tfuel - 1e-3 <= T <= Taw + 1e-3, ft.last_thermal_sol[2:end])
    end
end
