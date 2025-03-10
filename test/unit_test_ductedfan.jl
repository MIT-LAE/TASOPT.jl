@testset "Ducted fan" begin
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