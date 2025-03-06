
@testset "PEM fuel cell" begin

    u_LT = TASOPT.engine.PEMFC_inputs()

    u_LT.j = 1e4
    u_LT.T = 353.15
    u_LT.p_A = 3e5
    u_LT.p_C = 3e5
    u_LT.x_H2O_A = 0.1
    u_LT.x_H2O_C = 0.1
    u_LT.λ_H2 = 3.0
    u_LT.λ_O2 = 3.0
    u_LT.t_M = 100e-6
    u_LT.t_A = 250e-6
    u_LT.t_C = 250e-6
    u_LT.type = "LT-PEMFC"

    V_LT, α_LT = TASOPT.engine.LT_PEMFC_voltage(u_LT)

    @test V_LT ≈ 0.7031848337088604 
    @test α_LT ≈ 0.13026807266756726

    u_HT = TASOPT.engine.PEMFC_inputs()

    u_HT.j = 1e4
    u_HT.T = 453.15
    u_HT.p_A = 3e5
    u_HT.p_C = 3e5
    u_HT.x_H2O_A = 0.1
    u_HT.x_H2O_C = 0.1
    u_HT.λ_H2 = 3.0
    u_HT.λ_O2 = 3.0
    u_HT.t_M = 100e-6
    u_HT.t_A = 250e-6
    u_HT.t_C = 250e-6
    u_HT.type = "HT-PEMFC"

    V_HT= TASOPT.engine.HT_PEMFC_voltage(u_HT)

    @test V_HT ≈ 0.7511805292248419

    Pdes = 1e6
    Vdes = 200

    n_cells, A_cell, Q = TASOPT.engine.PEMsize(Pdes, Vdes, u_HT)

    @test n_cells ≈ 266.2475825969344
    @test A_cell ≈ 0.5
    @test Q ≈  676426.8383800443

    mfuel, V_stack, Q_od, j, α = TASOPT.engine.PEMoper(5e5, n_cells, A_cell, u_HT)

    @test mfuel ≈ 0.0058233133683716085 
    @test V_stack ≈ 238.82772890173237
    @test Q_od ≈ 201939.78148568503

    W_FC = TASOPT.engine.PEMstackweight(9.81, u_HT, n_cells, A_cell, 4.0)

    @test W_FC ≈ 7052.099720245001

    P2A_maxP, j_maxP = TASOPT.engine.find_maximum_PEMFC_power(u_LT)
    @test P2A_maxP ≈ 7614.221005014089
    @test j_maxP ≈ 13237.3046875

    W_FC, j = TASOPT.engine.PEM_weight_from_specific_power(u_LT, 1e7, 5e3, 2000)
    @test W_FC ≈ 74695.50805918823
    @test j ≈ j_maxP

end