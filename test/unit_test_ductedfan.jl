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

    out_size_check = (317.5455546477681, 0.45230025508427396, 3.175455546477681e6, 93.38924832690934, 245.96522147425617, -84555.91953868531, 34660.08023796853, 1004.6979001778267, 287.482334, 245.96522147425617, -84555.91953868531, 34660.08023796853, 1004.6979001778267, 287.482334, 245.96522147425617, -84555.91953868531, 34660.08023796853, 1004.6979001778267, 287.482334, 279.79114539158456, -50553.55199245569, 51990.120356952786, 1005.7466737051528, 287.482334, 279.79114539158456, -50553.55199245569, 51990.120356952786, 1005.7466737051528, 287.482334, 236.74253162629373, 229.40455573526395, 182.38072225868564, 27167.2752649851, 1004.3726660891869, 287.482334, 1.2430398905788191, 233.0847060242355, 306.3925547571069, 27455.55018043698, 1004.434294760138, 287.482334, 0.7438981200482442, 220.8920240309927, 344.04312976680137, 22756.739147503034, 1004.2508163054055, 287.482334, 0.7574695088798918, 0.9, 0.8941027193954981, true)

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

    out_opr_check = (139.5015648251009, 0.0, 57347.02696725949, 7.999999999999959e6, 227.10036959264153, 1.4349121928160546, 227.10036959264156, 1.0080619789152259, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 288.2, -42093.61205685147, 101320.0, 1006.0696710363969, 287.482334, 323.18144066294997, -6866.897039390857, 145385.30337612267, 1007.9050575955833, 287.482334, 323.18144066294997, -6866.897039390857, 145385.30337612267, 1007.9050575955833, 287.482334, 0.0, 273.89073896956006, 169.66987804142695, 84784.51444903799, 1005.5355852443635, 287.482334, 0.5109637742820122, 291.522920362439, 252.51842200928607, 101320.0, 1006.2071855389283, 287.482334, 0.7372058221722848, 291.522920362439, 252.51842200928607, 101320.0, 1006.2071855389283, 287.482334, 0.7372058221722849, 0.7438981200482441, 0.9, 0.8947714887329952)

    for (i,item) in enumerate(out_opr_takeoff) 
        @test item ≈ out_opr_check[i]
    end
end