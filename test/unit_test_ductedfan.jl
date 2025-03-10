@testset "Ducted fan" begin
    @testset "Ducetd fan models" begin
        #__ Test ducted fan sizing function __
        gee = TASOPT.gee
        M0 = 0.8
        h = 11.0 #km
        T0,p0,ρ0,a0,μ0 = TASOPT.atmos(h)
        M2 = 0.6
        Kinl = 0
        iBLIc = 0
        Phiinl = 0
        pifD = 1.5
        pid = 1.0
        pifn = 1.0
        epf0 = 0.9
        epfK = 0
        Δh_radiator = 0
        Δp_radiator = 0
        Fe = 1e4
        pifK = 0

        out_size = TASOPT.ductedfansize!(gee, M0, T0, p0, a0, M2,
            Fe, Phiinl, Kinl, iBLIc,
            pifD,
            pid, pifn, 
            epf0,
            pifK, epfK,
            Δh_radiator,
            Δp_radiator
            )

        TSEC, Fsp, Pfan, mfan,
        Tt0, ht0, pt0, cpt0, Rt0,
        Tt18, ht18, pt18, cpt18, Rt18,
        Tt2, ht2, pt2, cpt2, Rt2,
        Tt21, ht21, pt21, cpt21, Rt21,
        Tt7, ht7, pt7, cpt7, Rt7,
        u0,
        T2, u2, p2, cp2, R2, A2,
        T7, u7, p7, cp7, R7, A7,
        T8, u8, p8, cp8, R8, A8,
        epf,
        etaf,
        Lconv = out_size

        out_size_check = (323.25676181808547, 0.454271505100875, 3.2325676181808547e6, 92.98399826114989, 245.96522147425617, -84555.91953868531, 34660.08023796853, 1004.6979001778267, 287.482334, 245.96522147425617, -84555.91953868531, 34660.08023796853, 1004.6979001778267, 287.482334, 245.96522147425617, -84555.91953868531, 34660.08023796853, 1004.6979001778267, 287.482334, 280.54909072731016, -49791.14625907178, 51990.120356952786, 1005.7745681809906, 287.482334, 280.54909072731016, -49791.14625907178, 51990.120356952786, 1005.7745681809906, 287.482334, 236.74253162629373, 229.40455573526395, 182.38072225868564, 27167.2752649851, 1004.3726660891869, 287.482334, 1.237645886382155, 233.717293250894, 306.8073604765317, 27455.706502410256, 1004.4454671849026, 287.482334, 0.7416719063566514, 221.4913223688968, 344.50926624605864, 22756.739147503034, 1004.2585007339301, 287.482334, 0.755205520308143, 0.8814814814814815, 0.8744946801306193, true)
        
        for (i,item) in enumerate(out_size) 
            @test item ≈ out_size_check[i]
        end

        mbfD = mfan * sqrt(Tt2 / TASOPT.Tref) / (pt2 / TASOPT.pref)
        Nf = 1.0 #Arbitrarily set to 1 as only ratios matter
        Nbf = Nf / sqrt(Tt2 / TASOPT.Tref)
        NbfD = Nbf

        #__ Test ducted fan operation function __
        #First check that it provides the desired values at the design point
        out_opr_des  = TASOPT.ductedfanoper!(M0, T0, p0, a0, Tref, pref,
                    Phiinl, Kinl, iBLIc,
                    pid, pifn, 
                    pifD, 
                    mbfD,
                    NbfD, 
                    A2, A7,
                    epf0,
                    pifK, epfK,
                    Fe, 0,
                    M2, pifD, 0, 
                    Δh_radiator, Δp_radiator,
                    false)

        Fe2 = out_opr_des[3]
        Pf2 = out_opr_des[4]
        mfan2 = out_opr_des[5]
        pif2 = out_opr_des[6]

        #Test that the operation code returns the design variables
        @test Fe2 ≈ Fe
        @test Pf2 ≈ Pfan
        @test mfan2 ≈ mfan
        @test pif2 ≈ pifD

        #Now check conditions at takeoff 
        T0,p0,ρ0,a0,μ0 = TASOPT.atmos(0.0)
        M0 = 0
        Pf_takeoff = 8e6

        out_opr_takeoff  = TASOPT.ductedfanoper!(M0, T0, p0, a0, Tref, pref,
                    Phiinl, Kinl, iBLIc,
                    pid, pifn, 
                    pifD, 
                    mbfD,
                    NbfD, 
                    A2, A7,
                    epf0,
                    pifK, epfK,
                    0, Pf_takeoff,
                    M2, pifD, 0, 
                    Δh_radiator, Δp_radiator,
                    true)
        
        out_opr_check = (140.9379118321713, 0.0, 56762.58358028376, 8.000000000004287e6, 225.43165522123786, 1.4313614761026359, 225.4316552212379, 1.0048896010439536, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 323.44014553944726, -6606.1380720949155, 145025.54475871907, 1007.9214804675017, 287.482334, 323.44014553944726, -6606.1380720949155, 145025.54475871907, 1007.9214804675017, 287.482334, 0.0, 274.0073335407209, 168.97746825929116, 84910.82231481027, 1005.5396678345012, 287.482334, 0.5087707039345392, 291.9632588129839, 251.79508851397625, 101320.0, 1006.2259072821828, 287.482334, 0.7345423025364717, 291.9632588129839, 251.79508851397625, 101320.0, 1006.2259072821828, 287.482334, 0.7345423025364715, 0.7416719063568534, 0.8876187944374615, 0.8817854880614415)

        for (i,item) in enumerate(out_opr_takeoff) 
            @test item ≈ out_opr_check[i]
        end
    end

    @testset "Ducted fan with fuel cell" begin
        ac = TASOPT.load_default_model()
        pare = ac.pare
        parg = ac.parg

        #First, test the fuel cell + ducted fan design
        #Create a ducted fan with fuel cell model
        modelname = "fuel_cell_with_ducted_fan"
        engineweightname = "nasa"

        enginecalc! = TASOPT.engine.calculate_fuel_cell_with_ducted_fan!
        engineweight! = TASOPT.engine.fuel_cell_with_ducted_fan_weight!
        enginemodel = TASOPT.engine.FuelCellDuctedFan(modelname, enginecalc!, engineweightname, engineweight!)
        pare[iePfanmax,:,:] .= 20e6

        fcdata = TASOPT.engine.FuelCellDuctedFanData(2)

        fcdata.type = "HT-PEMFC"
        fcdata.current_density[iprotate,:] .= 1e4
        fcdata.FC_temperature .= 453.15
        fcdata.FC_pressure .= 3e5
        fcdata.water_concentration_anode .= 0.1
        fcdata.water_concentration_cathode .= 0.1
        fcdata.λ_H2 .= 3.0
        fcdata.λ_O2 .= 3.0
        fcdata.thickness_membrane = 100e-6
        fcdata.thickness_anode  = 250e-6
        fcdata.thickness_cathode  = 250e-6
        fcdata.design_voltage = 200.0
        pare[ieRadiatorepsilon,:,:] .= 0.7
        pare[ieRadiatorMp,:,:] .= 0.12
        pare[ieDi,:,:] .= 0.4

        para[iaROCdes, ipclimb1:ipclimbn,:] .= 500 * ft_to_m / 60
        engdata = fcdata

        engine = TASOPT.engine.Engine(enginemodel, engdata, Vector{TASOPT.engine.HX_struct}())

        ac.engine = engine

        #Prepare the pare object
        pare[ieRadiatorCoolantT,:,:] = engine.data.FC_temperature[:,:]
        pare[ieRadiatorCoolantP,:,:] = engine.data.FC_pressure[:,:]
        pare[ieRadiatorHeat,:,:] = engine.data.FC_heat[:,:]

        pare[ieFe,ipcruise1,1] = 16981.808185580507
        parg[igneng] = 2
        ac.wing.layout.S = 81.25043040696103
    
        pare[ieFsp,ipcruise1,1] = 0.5268888878557653
        pare[iepif,ipcruise1,1] = 1.685
        pare[iepid,ipcruise1,1] = 0.998
        pare[iepifn,ipcruise1,1] = 0.98
        
        pare[ieepolf,ipcruise1,1] = 0.8948
        
        pare[iepifK,ipcruise1,1] = 1.685
        pare[ieepfK,ipcruise1,1] = -0.077
        pare[ieM2,ipcruise1,1] = 0.6
        pare[ieM0,ipcruise1,1] = 0.8
        pare[ieTt0,ipcruise1,1] = 247.50937622665612
        pare[ieht0,ipcruise1,1] = -83004.69052610746
        pare[iept0,ipcruise1,1] = 36436.066979351635
        pare[iecpt0,ipcruise1,1] = 1004.734911318952
        pare[ieRt0,ipcruise1,1] = 287.482334
        pare[iep0,ipcruise1,1] = 23922.608843328788
        pare[iea0,ipcruise1,1] = 296.8557888469756
        pare[ierho0,ipcruise1,1] = 0.3800541947033053
        pare[iemu0,ipcruise1,1] = 1.4294279408521106e-5
        pare[ieT0,ipcruise1,1] = 219.43067572699252
        pare[ieu0,ipcruise1,1] = 237.4846310775805

        engine.enginecalc!(ac, "design", 1, ipcruise1, true)

        @test pare[iemfuel,ipcruise1,1] ≈ 0.0004947385392474482
        @test pare[iePfan,ipcruise1,1] ≈ 5.972673833105838e6
        @test ac.engine.data.number_cells ≈ 265.5192500533146
        @test ac.engine.data.area_cell ≈ 10.0
        @test ac.engine.data.FC_heat[ipcruise1,1] ≈ 1.9112072243161257e6

        #Next, test the off-design performance

        pare[ieFe,iprotate,1] = 85275.63659107959
        parg[igneng] = 2
        ac.wing.layout.S = 81.25043040696103
    
        pare[ieFsp,iprotate,1] = 3.25725288027002
        pare[iepif,iprotate,1] = 1.718876825147152
        pare[iepid,iprotate,1] = 0.998
        pare[iepifn,iprotate,1] = 0.98
        
        pare[ieepolf,iprotate,1] = 0.8948
        
        pare[iepifK,iprotate,1] = 1.685
        pare[ieepfK,iprotate,1] = -0.077
        pare[ieM2,iprotate,1] = 0.6
        pare[ieM0,iprotate,1] = 0.21721987853710417
        pare[ieTt0,iprotate,1] = 290.9134582811252
        pare[ieht0,iprotate,1] = -39363.02002852104
        pare[iept0,iprotate,1] = 104698.14961489625
        pare[iecpt0,iprotate,1] = 1006.18147266869
        pare[ieRt0,iprotate,1] = 287.482334
        pare[iep0,iprotate,1] = 101320.0
        pare[iea0,iprotate,1] = 340.2074661144284
        pare[ierho0,iprotate,1] = 1.2255627040761317
        pare[iemu0,iprotate,1] = 1.78e-5
        pare[ieT0,iprotate,1] = 288.2
        pare[ieu0,iprotate,1] = 73.89982446679213
        engine.enginecalc!(ac, "off_design", 1, iprotate, true)

        @test pare[iemfuel,iprotate,1] ≈ 0.002089436746114079
        @test pare[iePfan,iprotate,1] ≈ 2.0000000001539476e7
        @test ac.engine.data.FC_heat[iprotate,1] ≈ 1.3296113958830826e7
    end
end