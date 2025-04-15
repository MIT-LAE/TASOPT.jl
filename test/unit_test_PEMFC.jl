
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

    @test V_LT ≈ 0.7103339015901087
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

    @test V_HT ≈ 0.7532410548758377

    Pdes = 1e6
    Vdes = 200

    n_cells, A_cell, Q = TASOPT.engine.PEMsize(Pdes, Vdes, u_HT)

    @test n_cells ≈ 265.5192500533146
    @test A_cell ≈ 0.5
    @test Q ≈ 664805.6978342824

    V_stack, Q_od = TASOPT.engine.PEMoper(5e5, n_cells, A_cell, u_HT)

    @test V_stack ≈ 237.6830551160271
    @test Q_od ≈ 200430.95710865577

    W_FC = TASOPT.engine.PEMstackweight(9.81, u_HT, n_cells, A_cell, 4.0)

    @test W_FC ≈ 7032.808376162144
end