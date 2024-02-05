using TASOPT
using Test

@testset "PEM fuel cell" begin

    u_LT = TASOPT.engine.PEMFC_inputs(1e4, 353.15, 3e5, 3e5, 0.1, 0.1, 3.0, 3.0, 100e-6, 250e-6, 250e-6, "LT-PEMFC")

    V_LT, α_LT = TASOPT.engine.LT_PEMFC_voltage(u_LT)

    @test V_LT == 0.7102709836153386
    @test α_LT == 0.13026800553988963

    u_HT = TASOPT.engine.PEMFC_inputs(1e4, 453.15, 3e5, 3e5, 0.1, 0.1, 3.0, 3.0, 100e-6, 250e-6, 250e-6, "HT-PEMFC")
    V_HT= TASOPT.engine.HT_PEMFC_voltage(u_HT)

    @test V_HT == 0.7532410548758377

    Pdes = 1e6
    Vdes = 200

    n_cells, A_cell, Q = TASOPT.engine.PEMsize(Pdes, Vdes, u_HT)

    @test n_cells == 265.5192500533146
    @test A_cell == 0.5
    @test Q == 664805.6978342824

    V_stack, Q_od = TASOPT.engine.PEMoper(5e5, n_cells, A_cell, u_HT)

    @test V_stack == 237.6830551160271
    @test Q_od == 200430.95710865577

    W_FC = TASOPT.engine.PEMstackweight(9.81, u_HT, n_cells, A_cell, 4.0)

    @test W_FC == 7032.808376162144
end