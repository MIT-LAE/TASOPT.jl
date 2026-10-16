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

        @test inner_tank_mech.Wtank       ≈ 59963.9469836711
        @test inner_tank_mech.Winsul_sum  ≈ 38147.75085413718
        @test inner_tank_mech.Vfuel       ≈ 166.77327116515787
        @test inner_tank_mech.Shead_insul ≈ [16.421341718687554, 20.311098333705395, 24.62038329710079, 29.34797317883484]
        @test inner_tank_mech.Rtank_outer ≈ 1.874385890773477
        @test inner_tank_mech.l_tank      ≈ 15.771796698235905
        @test inner_tank_mech.l_cyl       ≈ 12.842691468445997
        
        outputs_thermal = TASOPT.CryoTank.tankWthermal(fuse, fuse_tank, z, TSL, Mair, xftank, ifuel)

        outputs_thermal_check = (1835.8610739583605,)

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

        @test inner_tank_vac.Wtank       ≈ 24003.269280331166
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
        @test outer_tank_vac.Wtank   ≈ 87480.86746096004
        @test outer_tank_vac.Wcyl    ≈ 31381.827647692957
        @test outer_tank_vac.Whead   ≈ 13786.613724211644
        @test outer_tank_vac.Wstiff  ≈ 20573.006232029245
        @test outer_tank_vac.Souter  ≈ 171.6223477616751
        @test outer_tank_vac.Shead   ≈ 26.90237940128641
        @test outer_tank_vac.Scyl    ≈ 117.81758895910225
        @test outer_tank_vac.t_cyl   ≈ 0.009560503271690743
        @test outer_tank_vac.t_head  ≈ 0.018394143361195915
        @test outer_tank_vac.l_outer ≈ 9.96235436442712

    end
    
    @testset "outer vessel uses outer_material instead of inner_material" begin
        ft = deepcopy(fuse_tank)
        ft.outer_material = StructuralAlloy("Ti-6242") # deliberately different from inner_material
        out_ti = TASOPT.CryoTank.size_outer_tank(fuse, ft, 7e4, 10.0, 1.0)
        ft.outer_material = StructuralAlloy("Al-2219-T87")
        out_al = TASOPT.CryoTank.size_outer_tank(fuse, ft, 7e4, 10.0, 1.0)
        @test out_ti.Wtank != out_al.Wtank # would be equal if inner_material were read
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

        # A composite inner wall takes its (higher) emissivity from the material
        Rvac_cf = TASOPT.CryoTank.vacuum_resistance(fuse_tank.Tfuel, Taw, 90.0, 100.0;
                        inner_material = TASOPT.materials.CompositeLaminate("CFRP-IM7-8552-PW-0090"))
        @test Rvac_cf < Rvac
        @test TASOPT.CryoTank.vacuum_resistance(fuse_tank.Tfuel, Taw, 90.0, 100.0;
                        inner_material = fuse_tank.inner_material) ≈ Rvac_check
    end
end

@testset "fuselage_tank composite fields (T3)" begin
    new_field_defaults = (wall_architecture = "monocoque", has_permeation_barrier = false,
        t_barrier = 0.0, f_joint_add = 0.0, design_strain_limit = 0.0,
        dome_layup_efficiency = 1.0, buckling_knockdown = 1.0, composite_min_stacks = 1,
        eta_composite = 1.0)

    @testset "new fields have inert defaults" begin
        ft = TASOPT.fuselage_tank()
        for (name, value) in pairs(new_field_defaults)
            @test getfield(ft, name) == value
        end
        @test ft.stiffener_material isa StructuralAlloy
        @test ft.stiffener_material.name == "Al-2219-T87"
    end

    @testset "wall materials accept CompositeLaminate or StructuralAlloy" begin
        ft = deepcopy(fuse_tank)
        cfrp = TASOPT.materials.CompositeLaminate("CFRP-IM7-8552-PW-0090")
        ft.inner_material = cfrp
        ft.outer_material = cfrp
        @test ft.inner_material isa TASOPT.materials.CompositeLaminate
        @test ft.outer_material isa TASOPT.materials.CompositeLaminate
        ft.inner_material = StructuralAlloy("Ti-6242")
        @test ft.inner_material isa StructuralAlloy
    end

    @testset "aluminium vacuum-tank sizing still matches the reference values" begin
        ft = deepcopy(fuse_tank)
        inner = TASOPT.CryoTank.size_inner_tank(fuse, ft, ft.t_insul)
        outer = TASOPT.CryoTank.size_outer_tank(fuse, ft, inner.Wtank + ft.Wfuelintank, inner.l_cyl, 14.025390625)
        @test inner.Wtank ≈ 24003.269280331166
        @test outer.Wtank ≈ 87480.86746096004
    end

    @testset "default model keeps an alloy inner wall" begin
        @test load_default_model().fuse_tank.inner_material isa StructuralAlloy
    end
end

@testset "Composite inner tank (T4)" begin
    CompositeLaminate = TASOPT.materials.CompositeLaminate
    cfrp0090 = CompositeLaminate("CFRP-IM7-8552-PW-0090")
    cfrpQI   = CompositeLaminate("CFRP-IM7-8552-PW-QI")
    Rout = 2.335 # outer radius of the inner vessel in the vacuum configuration
    wall(mat, ft) = TASOPT.CryoTank.composite_inner_wall(mat, ft, Rout)
    cfrp_tank(mat; kwargs...) = begin
        ft = deepcopy(fuse_tank)
        ft.inner_material = mat
        for (name, value) in kwargs
            setfield!(ft, name, value)
        end
        ft
    end
    size_W(ft) = TASOPT.CryoTank.size_inner_tank(fuse, ft, ft.t_insul).Wtank

    @testset "alloy wall ignores stiffener_material" begin
        ft = deepcopy(fuse_tank)
        ft.stiffener_material = StructuralAlloy("Ti-6242")
        @test size_W(ft) == size_W(fuse_tank)
    end

    @testset "composite_inner_wall" begin
        ft = cfrp_tank(cfrp0090)
        w = wall(cfrp0090, ft)
        n_stack0090 = cfrp0090.plies_per_stack # 5 plies in the stack [F0/F0/F45/F0/F0]
        n_stackQI   = cfrpQI.plies_per_stack   # 4 plies in the stack [F0/F45/F45/F0]

        @test w.n_cyl isa Int && w.n_head isa Int
        @test w.t_cyl  == w.n_cyl  * cfrp0090.t_ply
        @test w.t_head == w.n_head * cfrp0090.t_ply
        # the wall is built from whole stacks
        @test w.n_cyl % n_stack0090 == 0 && w.n_head % n_stack0090 == 0
        wqi = wall(cfrpQI, cfrp_tank(cfrpQI))
        @test wqi.n_cyl % n_stackQI == 0 && wqi.n_head % n_stackQI == 0
        @test wqi.t_cyl == wqi.n_cyl * cfrpQI.t_ply
        @test w.mode == :hoop
        @test !w.floor_governs

        # hand calculation: hoop load Δp R over E2 * strain allowable
        epsa = TASOPT.CryoTank.composite_strain_allowable(cfrp0090, ft) * cfrp0090.knockdown_cryo * cfrp0090.knockdown_process
        t_hoop = ft.pvent * Rout / (cfrp0090.E2 * epsa * ft.eta_composite)
        @test w.n_cyl == n_stack0090 * max(1, ceil(Int, t_hoop / (n_stack0090 * cfrp0090.t_ply) - 1e-9))

        # floor at low pressure: one stack by default, or the user-set minimum number of stacks
        wlow = wall(cfrp0090, cfrp_tank(cfrp0090; pvent = 1e3))
        @test wlow.floor_governs
        @test wlow.n_cyl == n_stack0090
        @test wall(cfrp0090, cfrp_tank(cfrp0090; pvent = 1e3, composite_min_stacks = 3)).n_cyl == 15
        @test wall(cfrpQI, cfrp_tank(cfrpQI; pvent = 1e3, composite_min_stacks = 3)).n_cyl == 12

        # at large Δp the ply count is large, so rounding (steps of 5 plies) does not hide the scaling;
        # the tolerance is widened to atol = 0.10 to cover that remaining step
        big = cfrp_tank(cfrp0090; pvent = 2e7)
        n0 = wall(cfrp0090, big)
        @test n0.n_head / n0.n_cyl ≈ 1.0 atol = 0.10 # sphere-like head (AR=2, K=1) sees the same load as the barrel hoop (the Δp terms of Barron's Eqs. 7.1/7.2 still differ by ~5 % at 2e7 Pa)
        @test wall(cfrp0090, cfrp_tank(cfrp0090; pvent = 2e7, dome_layup_efficiency = 0.9)).n_head / n0.n_head ≈ 1 / 0.9 atol = 0.10
        # a design strain 20 % below the material allowable needs ~1/0.8 = 1.25 times the plies
        @test wall(cfrp0090, cfrp_tank(cfrp0090; pvent = 2e7, design_strain_limit = 0.8 * TASOPT.CryoTank.composite_strain_allowable(cfrp0090, fuse_tank))).n_cyl / n0.n_cyl ≈ 1.25 atol = 0.10
        @test wall(cfrp0090, cfrp_tank(cfrp0090; dome_layup_efficiency = 0.9)).n_head >= w.n_head
    end

    @testset "full sizing" begin
        al   = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
        c0   = TASOPT.CryoTank.size_inner_tank(fuse,
                   cfrp_tank(cfrp0090;
                       design_strain_limit = 0.0025,    # ASSUMPTION - needs to be re-evaluated; 2500 microstrain chosen for a multi-use aircraft tank (one or more flights/day for 10-20 years, i.e. thousands of pressure cycles); the database value 5000 microstrain is the NASA CCTD launch-vehicle value
                       f_joint_add = 0.0,               # deliberately zero: no substantiated value exists (transparent zero rather than an unsourced number); the joint mass term is covered by the 'mass add-ons' testset
                       has_permeation_barrier = true,
                       t_barrier = 0*0.38e-3),            # ASSUMPTION - needs to be re-evaluated (m)
                   fuse_tank.t_insul)

        @test isfinite(c0.Wtank) && c0.Wtank > 0
        @test c0.Vfuel ≈ al.Vfuel
        @test c0.Rtank_outer ≈ al.Rtank_outer
        # Plan sanity upper bound is 60 % (above that something is probably missing, plan section 2.4.2). CCTD's ~30 %
        # is for the whole tank against Al-Li, whereas the baseline here is the plain UTS/4 aluminium rule.
        # Expected value with these inputs is ~52.6 % (15-ply barrel, three stacks of five).
        @test 0.10 < 1 - c0.Wtank / al.Wtank < 0.60

        # In the continuous limit (very large Δp, negligible ply quantisation) the 0090 wall is thinner than QI.
        # At the default pressure stack rounding (5 vs 4 plies) can reverse this (plan OQ3).
        @test wall(cfrp0090, cfrp_tank(cfrp0090; pvent = 2e7)).t_cyl < wall(cfrpQI, cfrp_tank(cfrpQI; pvent = 2e7)).t_cyl
    end

    @testset "mass add-ons" begin
        W0 = size_W(cfrp_tank(cfrp0090))
        # barrier thickness is inert unless the barrier is switched on
        @test size_W(cfrp_tank(cfrp0090; t_barrier = 0.4e-3)) == W0
        dWb1 = size_W(cfrp_tank(cfrp0090; has_permeation_barrier = true, t_barrier = 0.4e-3)) - W0
        dWb2 = size_W(cfrp_tank(cfrp0090; has_permeation_barrier = true, t_barrier = 0.8e-3)) - W0
        @test dWb1 > 0
        @test dWb2 ≈ 2 * dWb1 rtol = 1e-6

        dWj1 = size_W(cfrp_tank(cfrp0090; f_joint_add = 0.05)) - W0
        dWj2 = size_W(cfrp_tank(cfrp0090; f_joint_add = 0.10)) - W0
        @test dWj1 > 0
        @test dWj2 ≈ 2 * dWj1 rtol = 1e-6
    end

    @testset "ring material of a composite-wall tank" begin
        W_al = size_W(cfrp_tank(cfrp0090))
        W_ti = size_W(cfrp_tank(cfrp0090; stiffener_material = StructuralAlloy("Ti-6242")))
        @test isfinite(W_ti)
        @test W_ti != W_al
    end
end

@testset "Composite outer tank (T5)" begin
    CompositeLaminate = TASOPT.materials.CompositeLaminate
    cfrp0090 = CompositeLaminate("CFRP-IM7-8552-PW-0090")
    cfrpVDA  = CompositeLaminate("CFRP-IM7-8552-PW-0090-VDA")
    cfrpQI   = CompositeLaminate("CFRP-IM7-8552-PW-QI")
    Winner, l_cyl, Ninterm = 7e4, 10.0, 1.0 # inner-vessel load (N), cylinder length (m), intermediate rings
    outer_tank(mat; kwargs...) = begin
        ft = deepcopy(fuse_tank)
        ft.outer_material = mat
        for (name, value) in kwargs
            setfield!(ft, name, value)
        end
        ft
    end
    size_outer(ft; l = l_cyl, N = Ninterm) = TASOPT.CryoTank.size_outer_tank(fuse, ft, Winner, l, N)

    # Independent evaluation of the continuous required thicknesses (Barron 7.11 and the head formula)
    # for a wall with modulus E and Poisson ratio ν
    function required_thickness(ft, E, ν; l = l_cyl, N = Ninterm)
        Do = 2 * (fuse.layout.radius - ft.clearance_fuse)
        L_Do = l / (2 + N - 1) / Do
        pc = 4 * TASOPT.pref
        residual(t_D) = 2.42 * E * t_D^(5/2) / ((1 - ν^2)^(3/4) * (L_Do - 0.45 * sqrt(t_D))) - pc
        t_cyl = TASOPT.CryoTank.find_zero(residual, 1e-3) * Do
        K1 = TASOPT.CryoTank.find_K1_head(ft.ARtank)
        t_head = K1 * Do * sqrt(pc * sqrt(3 * (1 - ν^2)) / (0.5 * E))
        return t_cyl, t_head
    end
    t_stack(mat) = mat.plies_per_stack * mat.t_ply
    composite_required(ft; kwargs...) = required_thickness(ft,
        ft.outer_material.E_eq * ft.buckling_knockdown, ft.outer_material.ν_eq; kwargs...)

    @testset "alloy outer wall is unchanged" begin
        al = fuse_tank.outer_material
        out = size_outer(fuse_tank)
        # the independent requirement helper reproduces the existing Al-2219 outer wall
        t_cyl_req, t_head_req = required_thickness(fuse_tank, al.E, al.ν)
        @test out.t_cyl  ≈ t_cyl_req  rtol = 1e-6
        @test out.t_head ≈ t_head_req rtol = 1e-6
        # no facing on an alloy wall; the knockdown acts on composite walls only
        @test out.Wfacing == 0.0
        out_kd = size_outer(outer_tank(al; buckling_knockdown = 0.65))
        @test out_kd.Wtank == out.Wtank
        @test out_kd.t_cyl == out.t_cyl
        @test out_kd.t_head == out.t_head
        @test out_kd.Wfacing == 0.0
    end

    @testset "composite wall is built from whole stacks" begin
        ft = outer_tank(cfrp0090)
        out = size_outer(ft)
        t_cyl_req, t_head_req = composite_required(ft)
        ts = t_stack(cfrp0090)
        for (t, t_req) in ((out.t_cyl, t_cyl_req), (out.t_head, t_head_req))
            n_stacks = round(Int, t / ts)
            @test t ≈ n_stacks * ts rtol = 1e-12  # exact integer multiple of plies_per_stack * t_ply
            @test n_stacks >= ft.composite_min_stacks
            @test t >= t_req * (1 - 1e-9)
            # rounded up by less than one stack, unless the minimum-stack floor governs
            @test (t - t_req < ts) || (n_stacks == ft.composite_min_stacks)
        end
        # mass uses the cured laminate density
        @test out.Wcyl  ≈ out.Scyl  * out.t_cyl  * cfrp0090.ρ * TASOPT.gee rtol = 1e-9
        @test out.Whead ≈ out.Shead * out.t_head * cfrp0090.ρ * TASOPT.gee rtol = 1e-9
    end

    @testset "buckling knockdown and minimum stacks" begin
        out1  = size_outer(outer_tank(cfrp0090))
        out65 = size_outer(outer_tank(cfrp0090; buckling_knockdown = 0.65)) # ASSUMPTION - 0.65 is an illustrative knockdown
        @test out65.t_cyl  >= out1.t_cyl
        @test out65.t_head >= out1.t_head
        @test out65.Wtank  >= out1.Wtank
        ft65 = outer_tank(cfrp0090; buckling_knockdown = 0.65)
        t_cyl_req, t_head_req = composite_required(ft65)
        @test out65.t_cyl >= t_cyl_req * (1 - 1e-9)
        @test out65.t_head >= t_head_req * (1 - 1e-9)

        # a large minimum number of stacks forces the floor on both cylinder and head
        # ASSUMPTION - 40 stacks (about 40 mm) exceed the buckling requirement (of order 10 mm, similar to the Al-2219 jacket)
        out40 = size_outer(outer_tank(cfrp0090; composite_min_stacks = 40))
        @test out40.t_cyl  ≈ 40 * t_stack(cfrp0090) rtol = 1e-12
        @test out40.t_head ≈ 40 * t_stack(cfrp0090) rtol = 1e-12
        @test out40.Wtank > out1.Wtank
    end

    @testset "quasi-isotropic outer wall and comparison with aluminium" begin
        ftqi = outer_tank(cfrpQI)
        outqi = size_outer(ftqi)
        tsqi = t_stack(cfrpQI)
        @test outqi.t_cyl  ≈ round(Int, outqi.t_cyl  / tsqi) * tsqi rtol = 1e-12
        @test outqi.t_head ≈ round(Int, outqi.t_head / tsqi) * tsqi rtol = 1e-12
        t_cyl_req, t_head_req = composite_required(ftqi)
        @test outqi.t_cyl >= t_cyl_req * (1 - 1e-9)
        @test outqi.t_head >= t_head_req * (1 - 1e-9)

        # ASSUMPTION - with knockdown 1.0 the 0/90 composite jacket is lighter than the Al-2219 jacket
        @test size_outer(outer_tank(cfrp0090)).Wtank < size_outer(fuse_tank).Wtank
    end

    @testset "ring material of a composite outer wall" begin
        out_al = size_outer(outer_tank(cfrp0090))
        ft_ti = outer_tank(cfrp0090; stiffener_material = StructuralAlloy("Ti-6242"))
        out_ti = size_outer(ft_ti)
        @test out_ti.Wstiff != out_al.Wstiff
        @test out_ti.t_cyl == out_al.t_cyl # the ring alloy does not change the wall
        @test out_ti.t_head == out_al.t_head
        # rings use the stiffener material, not the laminate: the default (Al-2219) ring is much heavier per volume than CFRP would give
        ref_ring = TASOPT.CryoTank.size_outer_tank(fuse, outer_tank(StructuralAlloy("Al-2219-T87")), Winner, l_cyl, Ninterm)
        @test out_al.Wstiff ≈ ref_ring.Wstiff rtol = 1e-9
    end

    @testset "low-emissivity facing mass" begin
        outv = size_outer(outer_tank(cfrpVDA))
        @test cfrpVDA.facing_areal_mass > 0
        @test outv.Wfacing ≈ cfrpVDA.facing_areal_mass * outv.Souter * TASOPT.gee rtol = 1e-12
        @test outv.Wfacing > 0
        # included in the total weight: Wtank = (Wcyl + 2 Whead + Wstiff) (1 + ftankadd) + Wfacing
        @test outv.Wtank ≈ (outv.Wcyl + 2 * outv.Whead + outv.Wstiff) * (1 + fuse_tank.ftankadd) + outv.Wfacing rtol = 1e-12
        # the bare laminate carries no facing
        out0 = size_outer(outer_tank(cfrp0090))
        @test out0.Wfacing == 0.0
        @test out0.Wtank ≈ (out0.Wcyl + 2 * out0.Whead + out0.Wstiff) * (1 + fuse_tank.ftankadd) rtol = 1e-12
    end

    @testset "outer wall emissivity check" begin
        @test_logs (:warn,) TASOPT.CryoTank.check_outer_wall_emissivity(cfrp0090)
        @test_logs TASOPT.CryoTank.check_outer_wall_emissivity(cfrpVDA)
        @test_logs TASOPT.CryoTank.check_outer_wall_emissivity(StructuralAlloy("Al-2219-T87"))
    end

    @testset "optimize_outer_tank with a composite outer wall" begin
        ft = outer_tank(cfrp0090)
        Nopt = TASOPT.CryoTank.optimize_outer_tank(fuse, ft, Winner, l_cyl)
        @test isfinite(Nopt)
        @test 0.0 <= Nopt <= 50.0
        Wopt = size_outer(ft; N = Nopt).Wtank
        @test isfinite(Wopt) && Wopt > 0
    end
end

@testset "Outer-wall emissivity" begin
    CompositeLaminate = TASOPT.materials.CompositeLaminate
    Al = StructuralAlloy("Al-2219-T87")
    cfrp_bare = CompositeLaminate("CFRP-IM7-8552-PW-0090")
    cfrp_vda  = CompositeLaminate("CFRP-IM7-8552-PW-0090-VDA")
    Taw = 243.26644796151263
    Tcold = fuse_tank.Tfuel
    Rvac_al = 0.3565313051181731  # reference value, alloy inner and outer walls (eps = 0.04)

    # Independent re-implementation of the vacuum-gap model with separate inner/outer emissivities
    function rvac_reference(Tc, Th, Si, So, εi, εo)
        σ = TASOPT.σ_SB
        ai, ao, pv, Rgas, γ = 1.0, 0.9, 1e-2, 287.05, 1.4
        Fe = 1 / (1/εi + Si/So * (1/εo - 1))
        R_rad = 1 / (σ * Fe * (Tc^2 + Th^2) * (Tc + Th) * Si)
        Fa = 1 / (1/ai + Si/So * (1/ao - 1))
        G = (γ + 1)/(γ - 1) * sqrt(Rgas / (8*pi*Th)) * Fa
        R_conv = 1 / (G * pv * Si)
        return R_conv * R_rad / (R_conv + R_rad)
    end

    @testset "vacuum_resistance outer_material keyword" begin
        # ASSUMPTION: explicit alloy or nothing outer material keeps the eps = 0.04 path bit-identical
        @test TASOPT.CryoTank.vacuum_resistance(Tcold, Taw, 90.0, 100.0; inner_material = Al, outer_material = Al) == Rvac_al
        @test TASOPT.CryoTank.vacuum_resistance(Tcold, Taw, 90.0, 100.0; inner_material = Al, outer_material = nothing) == Rvac_al
        # ASSUMPTION: outer_material = nothing is treated like an alloy (eps = 0.04) for a composite inner wall too
        Rcf_nothing = TASOPT.CryoTank.vacuum_resistance(Tcold, Taw, 90.0, 100.0; inner_material = cfrp_bare, outer_material = nothing)
        Rcf_al      = TASOPT.CryoTank.vacuum_resistance(Tcold, Taw, 90.0, 100.0; inner_material = cfrp_bare, outer_material = Al)
        @test Rcf_nothing == Rcf_al
    end

    @testset "outer emissivity ordering" begin
        Rbare = TASOPT.CryoTank.vacuum_resistance(Tcold, Taw, 90.0, 100.0; inner_material = Al, outer_material = cfrp_bare)
        Rvda  = TASOPT.CryoTank.vacuum_resistance(Tcold, Taw, 90.0, 100.0; inner_material = Al, outer_material = cfrp_vda)
        # ASSUMPTION: low outer emissivity (less radiation) => larger vacuum resistance; bare CFRP (eps 0.88) leaks most
        @test Rbare < Rvda
        # VDA facing (eps 0.05) is close to the alloy (eps 0.04): slightly smaller resistance than the alloy outer wall
        @test Rvda < Rvac_al
        # ASSUMPTION: 10% tolerance is a loose ordering-level bound, not a calibrated value
        @test isapprox(Rvda, Rvac_al; rtol = 0.1)
    end

    @testset "emissivity factor matches independent formula" begin
        for (outer, eps_out) in ((cfrp_bare, cfrp_bare.ε_emis), (cfrp_vda, cfrp_vda.ε_emis))
            R = TASOPT.CryoTank.vacuum_resistance(Tcold, Taw, 90.0, 100.0; inner_material = Al, outer_material = outer)
            @test R ≈ rvac_reference(Tcold, Taw, 90.0, 100.0, 0.04, eps_out) rtol = 1e-10
        end
    end

    @testset "whole-tank heat leak with outer emissivity" begin
        TASOPT.tanksize!(ac, 1)  # the vacuum layer is set above; same sizing as "Vacuum insulation"
        Q_al = TASOPT.CryoTank.tankWthermal(fuse, fuse_tank, z, TSL, Mair, xftank, ifuel)
        ft_bare = deepcopy(fuse_tank)
        ft_bare.outer_material = cfrp_bare
        ft_vda = deepcopy(fuse_tank)
        ft_vda.outer_material = cfrp_vda
        Q_bare = TASOPT.CryoTank.tankWthermal(fuse, ft_bare, z, TSL, Mair, xftank, ifuel)
        Q_vda  = TASOPT.CryoTank.tankWthermal(fuse, ft_vda, z, TSL, Mair, xftank, ifuel)
        @test isfinite(Q_bare) && Q_bare > 0
        @test Q_bare > Q_vda  # bare CFRP outer vessel leaks more heat than the VDA-faced one
        # alloy outer wall: unchanged from the baseline
        ft_al = deepcopy(fuse_tank)
        ft_al.outer_material = Al
        @test TASOPT.CryoTank.tankWthermal(fuse, ft_al, z, TSL, Mair, xftank, ifuel) == Q_al
    end
end

@testset "Wall thermal resistance" begin
    CompositeLaminate = TASOPT.materials.CompositeLaminate
    cfrp0090 = CompositeLaminate("CFRP-IM7-8552-PW-0090")
    Al = StructuralAlloy("Al-2219-T87")

    @testset "wall_conductivity and wall_resistance" begin
        @test TASOPT.CryoTank.wall_conductivity(Al) == 121.0
        @test TASOPT.CryoTank.wall_conductivity(Al) == Al.k
        @test TASOPT.CryoTank.wall_conductivity(cfrp0090) == 0.8
        @test TASOPT.CryoTank.wall_conductivity(cfrp0090) == cfrp0090.k_tt
        # ASSUMPTION: a missing (k_tt = 0) CompositeLaminate conductivity falls back to Al-2219 (121 W/m/K)
        cfrp_nok = deepcopy(cfrp0090)
        cfrp_nok.k_tt = 0.0
        @test TASOPT.CryoTank.wall_conductivity(cfrp_nok) == 121.0
        @test TASOPT.CryoTank.wall_resistance(cfrp0090, 0.01, 50.0) ≈ 0.01 / (0.8 * 50.0)
        @test TASOPT.CryoTank.wall_resistance(Al, 4.0e-3, 90.0) ≈ 4.0e-3 / (121.0 * 90.0)
    end

    @testset "new fuselage_tank wall fields and InnerTankGeometry.t_cyl" begin
        ft0 = TASOPT.fuselage_tank()
        @test ft0.t_wall_inner == 0.0
        @test ft0.t_wall_outer == 0.0

        inner_al = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
        @test inner_al.t_cyl > 0.0

        ft_cf = deepcopy(fuse_tank)
        ft_cf.inner_material = cfrp0090
        inner_cf = TASOPT.CryoTank.size_inner_tank(fuse, ft_cf, ft_cf.t_insul)
        # ASSUMPTION: a composite inner wall is built from whole plies, so t_cyl is an integer number of plies
        @test inner_cf.t_cyl == TASOPT.CryoTank.composite_inner_wall(cfrp0090, ft_cf, inner_cf.Rtank_outer).t_cyl
        n_plies = inner_cf.t_cyl / cfrp0090.t_ply
        @test isapprox(n_plies, round(n_plies); atol = 1e-6)
    end

    @testset "tanksize! stores the sized wall thicknesses" begin
        TASOPT.tanksize!(ac, 1)  # vacuum configuration set at the top of this file
        inner = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
        @test fuse_tank.t_wall_inner > 0.0
        @test fuse_tank.t_wall_inner ≈ inner.t_cyl
        # outer vessel is sized on the scaled length, as in tanksize!
        Routertank = fuse.layout.radius - fuse_tank.clearance_fuse
        lcyl2 = inner.l_cyl * Routertank / inner.Rtank_outer
        outer = TASOPT.CryoTank.size_outer_tank(fuse, fuse_tank, inner.Wtank + fuse_tank.Wfuelintank, lcyl2, fuse_tank.Ninterm)
        @test fuse_tank.t_wall_outer > 0.0
        @test fuse_tank.t_wall_outer ≈ outer.t_cyl
    end

    @testset "wall resistance in the heat-leak model" begin
        TASOPT.tanksize!(ac, 1)
        heat(ft) = TASOPT.CryoTank.tankWthermal(fuse, ft, z, TSL, Mair, xftank, ifuel)
        # ASSUMPTION: t_wall = 0.0 means "not yet sized" and the thermal model uses 4 mm (nominal Al)
        ft_default = deepcopy(fuse_tank)
        ft_default.t_wall_inner = 0.0 # not yet known => 4 mm assumed
        ft_default.t_wall_outer = 0.0
        ft_4mm = deepcopy(fuse_tank)
        ft_4mm.t_wall_inner = 4.0e-3
        ft_4mm.t_wall_outer = 4.0e-3
        Q_default = heat(ft_default)
        @test Q_default == heat(ft_4mm)
        @test isfinite(Q_default) && Q_default > 0

        # ASSUMPTION: a thick CFRP wall (k_tt = 0.8) adds a real resistance, so the heat leak falls
        cfrp0090 = TASOPT.materials.CompositeLaminate("CFRP-IM7-8552-PW-0090")
        inner_thin = deepcopy(fuse_tank); inner_thin.inner_material = cfrp0090; inner_thin.t_wall_inner = 1e-5
        inner_thick = deepcopy(fuse_tank); inner_thick.inner_material = cfrp0090; inner_thick.t_wall_inner = 0.1
        Q_thin = heat(inner_thin)
        Q_thick = heat(inner_thick)
        @test isfinite(Q_thin) && Q_thin > 0
        @test isfinite(Q_thick) && Q_thick > 0
        @test Q_thick < Q_thin

        outer_thin = deepcopy(fuse_tank); outer_thin.outer_material = cfrp0090; outer_thin.t_wall_outer = 1e-5
        outer_thick = deepcopy(fuse_tank); outer_thick.outer_material = cfrp0090; outer_thick.t_wall_outer = 0.1
        @test heat(outer_thick) < heat(outer_thin)
    end

    @testset "outer-wall emissivity reaches the insulation sizing residual" begin
        # ASSUMPTION: the vacuum layer is the only segment, so iinsuldes = [1] (res_MLI_thick would index past t_insul otherwise)
        cfrp_bare = TASOPT.materials.CompositeLaminate("CFRP-IM7-8552-PW-0090")
        cfrp_vda  = TASOPT.materials.CompositeLaminate("CFRP-IM7-8552-PW-0090-VDA")
        TASOPT.tanksize!(ac, 1)
        function res_for(outer)
            ft = deepcopy(fuse_tank)
            ft.outer_material = outer
            ft.iinsuldes = [1]
            x = [0.0, fuse_tank.Tfuel + 1.0, fuse_tank.Tfuel + 100.0]  # [Δt, T_w, T_fuse] for one vacuum layer; arbitrary but fixed
            TASOPT.CryoTank.res_MLI_thick(x, fuse, ft, z, TSL, Mair, xftank, ifuel)
        end
        # ASSUMPTION: with the T7A gap-fix the outer wall emissivity enters the vacuum layer of the sizing residual
        @test res_for(cfrp_bare) != res_for(cfrp_vda)
    end
end

@testset "Fuel-dependent strain allowable" begin
    CompositeLaminate = TASOPT.materials.CompositeLaminate
    cfrp0090 = CompositeLaminate("CFRP-IM7-8552-PW-0090")
    cfrpQI   = CompositeLaminate("CFRP-IM7-8552-PW-QI")
    cfrpVDA  = CompositeLaminate("CFRP-IM7-8552-PW-0090-VDA")
    Rout = 2.335 # outer radius of the inner vessel in the vacuum configuration
    fuel_tank(mat, fuel; pvent = fuse_tank.pvent) = begin
        ft = deepcopy(fuse_tank)
        ft.inner_material = mat
        ft.fueltype = fuel
        ft.design_strain_limit = 0.0
        ft.pvent = pvent
        ft
    end
    strain(mat, fuel) = TASOPT.CryoTank.composite_strain_allowable(mat, fuel_tank(mat, fuel))

    @testset "fuel table values" begin
        @test strain(cfrp0090, "LH2")  ≈ 1.4e-3 rtol = 1e-12
        @test strain(cfrp0090, "CH4")  ≈ 1.8e-3 rtol = 1e-12
        @test strain(cfrp0090, "C2H4") ≈ 2.0e-3 rtol = 1e-12
        @test strain(cfrp0090, "C2H6") ≈ 2.1e-3 rtol = 1e-12
        # same design strain for all CFRP laminates
        @test strain(cfrpQI, "LH2")  ≈ 1.4e-3 rtol = 1e-12
        @test strain(cfrpVDA, "LH2") ≈ 1.4e-3 rtol = 1e-12
        # fuel matching is case-insensitive
        @test strain(cfrp0090, "lh2") ≈ 1.4e-3 rtol = 1e-12
    end

    @testset "design_strain_limit overrides the table; unknown fuel throws" begin
        ft = fuel_tank(cfrp0090, "LH2")
        ft.design_strain_limit = 0.0025
        @test TASOPT.CryoTank.composite_strain_allowable(cfrp0090, ft) == 0.0025
        @test_throws ArgumentError TASOPT.CryoTank.composite_strain_allowable(cfrp0090, fuel_tank(cfrp0090, "XYZ"))
        @test_throws ArgumentError TASOPT.CryoTank.composite_strain_allowable(cfrp0090, fuel_tank(cfrp0090, ""))
    end

    @testset "fuel sets the composite wall thickness" begin
        # continuous limit (pvent = 2e7) so that ply quantisation does not hide the scaling
        wall(mat, fuel) = TASOPT.CryoTank.composite_inner_wall(mat, fuel_tank(mat, fuel; pvent = 2e7), Rout)
        wLH2, wCH4, wC2H4, wC2H6 = (wall(cfrp0090, f) for f in ("LH2", "CH4", "C2H4", "C2H6"))
        # lower allowable strain needs a thicker wall
        @test wLH2.n_cyl >= wCH4.n_cyl >= wC2H4.n_cyl >= wC2H6.n_cyl
        @test wLH2.n_head >= wCH4.n_head >= wC2H4.n_head >= wC2H6.n_head
        # ratio of the continuous ply counts follows the ratio of the allowables (2.1 / 1.4)
        @test wLH2.n_cyl / wC2H6.n_cyl ≈ 2.1 / 1.4 atol = 0.10
    end
end

@testset "Sizing diagnostics" begin
    # Diagnostics the sizing functions return next to the wall thicknesses:
    # ply counts, wall thicknesses, governing mode per wall (cylinder and head). Alloy walls report :alloy.
    CompositeLaminate = TASOPT.materials.CompositeLaminate
    cfrp0090 = CompositeLaminate("CFRP-IM7-8552-PW-0090")
    Rout = 2.335 # outer radius of the inner vessel in the vacuum configuration
    Winner, l_cyl, Ninterm = 7e4, 10.0, 1.0 # inner-vessel load (N), cylinder length (m), intermediate rings
    cfrp_inner(; kwargs...) = begin
        ft = deepcopy(fuse_tank)
        ft.inner_material = cfrp0090
        for (name, value) in kwargs
            setfield!(ft, name, value)
        end
        ft
    end
    outer_tank(mat; kwargs...) = begin
        ft = deepcopy(fuse_tank)
        ft.outer_material = mat
        for (name, value) in kwargs
            setfield!(ft, name, value)
        end
        ft
    end
    size_outer(ft) = TASOPT.CryoTank.size_outer_tank(fuse, ft, Winner, l_cyl, Ninterm)

    @testset "composite inner wall modes" begin
        # high pressure (2e7 Pa, continuous limit): the barrel is set by hoop or axial strength, the head by strength
        big = TASOPT.CryoTank.composite_inner_wall(cfrp0090, cfrp_inner(pvent = 2e7), Rout)
        @test big.mode_cyl in (:hoop, :axial)
        @test big.mode_cyl == big.mode
        @test big.mode_head == :strength
        @test big.n_cyl ≈ big.t_cyl / cfrp0090.t_ply rtol = 1e-12
        @test big.n_cyl % cfrp0090.plies_per_stack == 0
        # default pressure: the floor does not govern, so the barrel mode is the strength mode
        w0 = TASOPT.CryoTank.composite_inner_wall(cfrp0090, cfrp_inner(), Rout)
        @test !w0.floor_governs
        @test w0.mode_cyl == w0.mode
        # low pressure: the stack floor (one stack by default) sets both walls
        low = TASOPT.CryoTank.composite_inner_wall(cfrp0090, cfrp_inner(pvent = 1e3), Rout)
        @test low.floor_governs
        @test low.mode_cyl == :min_stack
        @test low.mode_head == :min_stack
    end

    @testset "inner wall diagnostics in the tank geometry" begin
        # composite inner wall: ply counts and thicknesses match the wall, modes come from composite_inner_wall
        ft = cfrp_inner(pvent = 2e7)
        inner = TASOPT.CryoTank.size_inner_tank(fuse, ft, ft.t_insul)
        d = inner.wall_diag
        @test d isa TASOPT.CryoTank.WallDiagnostics
        @test d.n_cyl > 0 && d.n_head > 0
        @test d.t_cyl ≈ inner.t_cyl rtol = 1e-12
        @test d.n_cyl * cfrp0090.t_ply ≈ inner.t_cyl rtol = 1e-12
        @test d.mode_cyl in (:hoop, :axial)
        @test d.mode_head == :strength
        # default alloy inner wall: no plies, both modes :alloy, thickness equal to the geometry
        ial = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
        @test ial.wall_diag.n_cyl == 0 && ial.wall_diag.n_head == 0
        @test ial.wall_diag.mode_cyl == :alloy && ial.wall_diag.mode_head == :alloy
        @test ial.wall_diag.t_cyl == ial.t_cyl
        @test ial.wall_diag.t_head > 0 # InnerTankGeometry has no t_head field; thickness checked via wall_diag
    end

    @testset "outer wall diagnostics in the tank geometry" begin
        # buckling governs the composite jacket at the default one-stack floor
        out = size_outer(outer_tank(cfrp0090))
        d = out.wall_diag
        @test d isa TASOPT.CryoTank.WallDiagnostics
        @test d.mode_cyl == :buckling
        @test d.n_cyl * cfrp0090.t_ply ≈ out.t_cyl rtol = 1e-12
        @test d.t_cyl == out.t_cyl
        @test d.n_head * cfrp0090.t_ply ≈ out.t_head rtol = 1e-12
        # a large stack floor forces :min_stack on both walls
        out40 = size_outer(outer_tank(cfrp0090; composite_min_stacks = 40))
        @test out40.wall_diag.mode_cyl == :min_stack
        @test out40.wall_diag.mode_head == :min_stack
        # alloy outer wall: no plies, modes :alloy, thickness is the required thickness
        outal = size_outer(fuse_tank)
        @test outal.wall_diag.n_cyl == 0 && outal.wall_diag.n_head == 0
        @test outal.wall_diag.mode_cyl == :alloy && outal.wall_diag.mode_head == :alloy
        @test outal.wall_diag.t_cyl == outal.t_cyl
    end

    @testset "alloy sizing values unchanged by the diagnostics" begin
        # same reference values as 'aluminium vacuum-tank sizing still matches the reference values'
        ft = deepcopy(fuse_tank)
        inner = TASOPT.CryoTank.size_inner_tank(fuse, ft, ft.t_insul)
        outer = TASOPT.CryoTank.size_outer_tank(fuse, ft, inner.Wtank + ft.Wfuelintank, inner.l_cyl, 14.025390625)
        @test inner.Wtank ≈ 24003.269280331166
        @test outer.Wtank ≈ 87480.86746096004
    end
end
