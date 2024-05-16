#----------------------------------------------------------------
#     Heat exchanger sizing tests
#
#     Thos script determines the geometric parameters for a cross-flow heat 
#     exchanger in air. It assumes that air in flowing through the exchanger
#     and heat is absorbed by a coolant fluid.

#---------------------------------     
# Test functions individually
#---------------------------------

@testset "Heat Exchanging" begin
    @testset "HEX functions" begin
        #---------------------------------     
        # hxsize!()
        #---------------------------------
        HXgas = TASOPT.engine.HX_gas()
        HXgeom = TASOPT.engine.HX_tubular()

        HXgas.fluid_p = "air"
        HXgas.fluid_c = "h2"
        HXgas.mdot_p = 1144/60
        HXgas.mdot_c = 9.95/60
        HXgas.ε = 0.8
        HXgas.Tp_in = 778
        HXgas.Tc_in = 264
        HXgas.Mp_in  = 0.19
        HXgas.Mc_in = 0.0285
        HXgas.pp_in = 40e3
        HXgas.pc_in = 1515e3

        HXgas.alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        HXgas.igas_c = 40

        HXgeom.fconc = 1
        HXgeom.frecirc = 0
        HXgeom.D_i = 0.564
        HXgeom.l = 0.6084530646014857 #tube length
        HXgeom.n_stages = 4
        HXgeom.xt_D = 6
        HXgeom.xl_D = 1.25
        HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
        HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")
        HXgeom.Δpdes = 3e6

        TASOPT.engine.hxsize!(HXgas, HXgeom)

        size_out = [HXgas.Tp_out, HXgas.Tc_out, HXgas.Δp_p, HXgeom.N_t, HXgeom.n_passes, HXgeom.tD_o, HXgeom.A_cs]

        size_out_check = 
        [731.5888605437423, 665.8848846504773, 1386.3503589746103, 62.03322510460286, 8.04559242739893, 0.004760508726403918, 1.0189779296746375]

        for i = 1: length(size_out)
            @test size_out[i] ≈ size_out_check[i]
        end
        #---------------------------------     
        # hxweight()
        #---------------------------------
        gee = 9.81
        fouter = 1

        W = TASOPT.engine.hxweight(gee, HXgeom, fouter)

        W_check = 801.5192810553101

        @test W == W_check
        #---------------------------------     
        # hxoper!()
        #---------------------------------
        HXgas = TASOPT.engine.HX_gas()
        HXgeom = TASOPT.engine.HX_tubular()

        HXgas.fluid_p = "air"
        HXgas.fluid_c = "h2"
        HXgas.mdot_p = 2*1144/60
        HXgas.mdot_c = 2*9.95/60
        HXgas.ε = 0.8
        HXgas.Tp_in = 778
        HXgas.Tc_in = 264
        HXgas.Mp_in  = 0.19
        HXgas.Mc_in = 0.0285
        HXgas.pp_in = 3*40e3
        HXgas.pc_in = 3*1515e3

        HXgas.alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        HXgas.igas_c = 40

        HXgeom.fconc = 1
        HXgeom.frecirc = 0
        HXgeom.D_i = 0.564
        HXgeom.t = 0.03e-2 #m, wall thicknesss
        HXgeom.tD_o = 0.004760326082769499
        HXgeom.A_cs = 1.0189779296746375
        HXgeom.l = 0.6084530646014857 #tube length
        HXgeom.n_stages = 4
        HXgeom.n_passes = 8
        HXgeom.N_t = 62
        HXgeom.xt_D = 6
        HXgeom.xl_D = 1.25
        HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
        HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")

        TASOPT.engine.hxoper!(HXgas, HXgeom)

        oper_out = [HXgas.Tp_out, HXgas.Tc_out, HXgas.Δp_p, HXgas.ε]

        oper_out_check = [740.3160720471974, 591.0871550274949, 1657.0074667280992, 0.6501725929032108]

        for i = 1: length(oper_out)
            @test oper_out[i] ≈ oper_out_check[i]
        end

        #---------------------------------     
        # hxoptim!()
        #---------------------------------
        HXgas = TASOPT.engine.HX_gas()
        HXgeom = TASOPT.engine.HX_tubular()

        HXgas.fluid_p = "air"
        HXgas.fluid_c = "h2"
        HXgas.mdot_p = 1144/60
        HXgas.mdot_c = 9.95/60
        HXgas.ε = 0.8
        HXgas.Tp_in = 778
        HXgas.Tc_in = 264
        HXgas.Mp_in  = 0.19
        HXgas.pp_in = 40e3
        HXgas.pc_in = 1515e3

        HXgas.alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        HXgas.igas_c = 40

        HXgeom.fconc = 1
        HXgeom.frecirc = 0
        HXgeom.D_i = 0.564
        HXgeom.l = 0.6084530646014857 #tube length
        HXgeom.xl_D = 1
        HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
        HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")
        HXgeom.Δpdes = 3e6

        #Calculate starting point
        #First calculate minimum tube length
        _, _, _, _, cp_p_in, Rp = TASOPT.engine.gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)
        γ_p_in = cp_p_in / (cp_p_in - Rp)
        ρ_p_in = HXgas.pp_in / (Rp * HXgas.Tp_in)
        Vp_in = HXgas.Mp_in * sqrt(γ_p_in * Rp * HXgas.Tp_in)

        A_cs = HXgas.mdot_p / (ρ_p_in * Vp_in) #Cross-sectional area of freestream

        if HXgeom.fconc #Flow is concentric
            D_i = HXgeom.D_i
            D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

            lmin = (D_o - D_i) / 2 #minimum tube length
            linit = 1.1 * lmin
        else #square cross-section
            AR_min = 0.1 #Minimum aspect ratio
            lmin = sqrt(AR_min * A_cs)
            linit = sqrt(A_cs)
        end

        #Now set starting point
        initial_x = [3, 4, 4, linit] #Initial guess
        
        # Optimize heat exchanger design parameters
        TASOPT.engine.hxoptim!(HXgas, HXgeom, initial_x)

        optim_out = [HXgas.Mc_in, HXgeom.n_stages, HXgeom.xt_D, HXgeom.l]

        TASOPT.engine.hxsize!(HXgas, HXgeom)

        Iobj = HXgas.Pl_p + HXgas.Pl_c #Optimizer may choose slightly different points with similar objective function. Check I too
        I_check = 71956.07072651063

        @test Iobj ≈ I_check

        #---------------------------------     
        # hxoptim!() with recirculation and rectangular
        #---------------------------------
        HXgas = TASOPT.engine.HX_gas()
        HXgeom = TASOPT.engine.HX_tubular()

        HXgas.fluid_p = "air"
        HXgas.fluid_c = "h2"
        HXgas.mdot_p = 2 * 49.9/60
        HXgas.mdot_c = 9.95/60
        HXgas.ε = 0.825
        HXgas.Tp_in = 791
        HXgas.Tc_in = 20
        HXgas.Mp_in  = 0.01
        HXgas.pp_in = 1515e3
        HXgas.pc_in = 1515e3
        HXgas.recircT  = 200
        HXgas.h_lat  = 446e3 + 670e3

        HXgas.alpha_p = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127]
        HXgas.igas_c = 40

        HXgeom.fconc = 0
        HXgeom.frecirc = 1
        HXgeom.t = 0.03e-2 #m, wall thicknesss
        HXgeom.xl_D = 1
        HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W
        HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")
        HXgeom.Δpdes = 3e6

        #Calculate starting point
        #First calculate minimum tube length
        _, _, _, _, cp_p_in, Rh = TASOPT.engine.gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)
        γ_p_in = cp_p_in / (cp_p_in - Rh)
        ρ_p_in = HXgas.pp_in / (Rh * HXgas.Tp_in)
        Vp_in = HXgas.Mp_in * sqrt(γ_p_in * Rh * HXgas.Tp_in)

        A_cs = HXgas.mdot_p / (ρ_p_in * Vp_in) #Cross-sectional area of freestream

        if HXgeom.fconc #Flow is concentric
            D_i = HXgeom.D_i
            D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

            lmin = (D_o - D_i) / 2 #minimum tube length
            linit = 1.1 * lmin
        else #square cross-section
            AR_min = 0.1 #Minimum aspect ratio
            lmin = sqrt(AR_min * A_cs)
            linit = 1.1 * lmin
        end

        #Now set starting point
        initial_x = [2, 8, 2, linit] #Initial guess

        TASOPT.engine.hxoptim!(HXgas, HXgeom, initial_x)

        TASOPT.engine.hxsize!(HXgas, HXgeom)

        optimrec_out = [HXgas.Mc_in, HXgas.mdot_r, HXgeom.l, HXgeom.n_stages, HXgeom.xt_D]

        Iobj_rec = HXgas.Pl_p + HXgas.Pl_c #Optimizer may choose slightly different points with similar objective function. 

        I_check_rec = 3489.323299331222 

        @test Iobj_rec ≈ I_check_rec
    end

    @testset "Thermal and pressure models" begin
        #Colburn j-factor
        Re = 1e4
        j,Cf = TASOPT.engine.jcalc_pipe(Re)

        Cf_check = 0.0791 * Re^(-0.25)
        j_check = Cf_check/2
        @test Cf ≈ Cf_check
        @test j ≈ j_check

        #Nu past tubes
        Res = [20, 500, 1e4, 3e5]
        Pr = 0.7
        N_L = 10.0
        xt_D = 3.0
        xl_D = 1.0
        Nu_check = [2.970870635447898, 13.683013469507236, 86.34877055350766, 787.9216560234382]

        for (i,Re) in enumerate(Res)
            Nu = TASOPT.engine.Nu_calc_staggered_cyl(Re, Pr, N_L, xt_D, xl_D)
            @test Nu ≈ Nu_check[i]
        end

        #Δp past tubes
        Re = 25e3
        Acs = pi * (1.265^2 - 0.564^2)/4
        l = 0.6
        tD_o = 4.78e-3
        L = 0.19
        N_tubes_tot = 1984
        xt_D = 6.0
        xl_D = 1.25
        μ_μw = 1.3

        G = 1144/60 / (Acs - 62*tD_o*l) #From Brewer
        NFV = Acs * L - N_tubes_tot * pi * tD_o^2 * l / 4 #Net free volume
        Ah = N_tubes_tot * tD_o * pi * l
        Dv = 4 * NFV / Ah
        ρ = 40e3 / (287 * 750)
        
        Δp = TASOPT.engine.Δp_calc_staggered_cyl(Re, G, L, ρ, Dv, tD_o, xt_D, xl_D, μ_μw)
        Δp_check = 1323.2936826807222

        @test Δp ≈ Δp_check

        #Tube sizing
        HXgeom = TASOPT.engine.HX_tubular()
        HXgeom.material = TASOPT.StructuralAlloy("SS-304")

        Acc = 1.372e-5 * 4 * 62
        b = pi *  0.564
        n_stages = 4
        K = pi * b * n_stages / (4 * xt_D * Acc)
        HXgeom.Δpdes = 3e6
        TASOPT.engine.tubesize!(K, HXgeom)

        t_check = 300e-6
        tDo_check = 0.004792450000885403
        @test HXgeom.t ≈ t_check
        @test HXgeom.tD_o ≈ tDo_check

        K = 10 #different case with wall thickness above minimum
        t_check = 0.0014766145463702754
        tDo_check = 0.10582404248986974

        TASOPT.engine.tubesize!(K, HXgeom)
        @test HXgeom.t ≈ t_check
        @test HXgeom.tD_o ≈ tDo_check

    end

end