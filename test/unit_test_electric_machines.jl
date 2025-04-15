using TASOPT
using Test

@testset "Electric machines" begin
    #Motor
    mot = TASOPT.propsys.ElectricMachine.Motor()
    mot.N_pole_pairs = 10
    mot.J_max = 10e6
    mot.N_inverters = 10
    mot.windings.kpf = 0.35
    mot.U_max = 200.0
    mot.N_inverters = 1

    Pmech = 1e6
    shaft_speed = 10e3 #RPM
    TASOPT.propsys.ElectricMachine.size_PMSM!(mot, shaft_speed, Pmech)

    @test mot.mass ≈ 99.05953856907007 rtol = 1e-8
    @test mot.l ≈ 0.15181949029252456 rtol = 1e-8
    @test mot.P_input ≈ 1.033907198036092e6 rtol = 1e-8
    @test mot.phase_resistance ≈ 0.0009936903085147378 rtol = 1e-8

    Pmech = 2.5e5
    shaft_speed = 5e3
    TASOPT.propsys.ElectricMachine.operate_PMSM!(mot, shaft_speed, Pmech)
    @test mot.P_input ≈ 259230.26242272905 rtol = 1e-8

    #Generator
    gen = TASOPT.propsys.ElectricMachine.Generator()
    gen.N_pole_pairs = 10
    gen.J_max = 10e6
    gen.windings.kpf = 0.35
    gen.U_max = 200.0

    Pmech = 1e6
    shaft_speed = 10e3 #RPM
    TASOPT.propsys.ElectricMachine.size_PMSM!(gen, shaft_speed, Pmech)

    @test gen.mass ≈ 99.05953856907007 rtol = 1e-8
    @test gen.l ≈ 0.15181949029252456 rtol = 1e-8
    @test gen.P_output ≈ 966092.8019639081 rtol = 1e-8
    @test gen.phase_resistance ≈ 0.0009936903085147378 rtol = 1e-8

    Pmech = 2.5e5
    shaft_speed = 5e3
    TASOPT.propsys.ElectricMachine.operate_PMSM!(gen, shaft_speed, Pmech)
    @test gen.P_output ≈ 240769.73757727098 rtol = 1e-8

    #Inverter
    invtr = TASOPT.propsys.ElectricMachine.Inverter()
    P_oupt = 1e6
    f = 800.0
    TASOPT.propsys.ElectricMachine.size_inverter!(invtr,P_oupt, f)

    @test invtr.mass ≈ 52.63157894736842 rtol = 1e-8
    @test invtr.P_input ≈ 1.0090817356205853e6 rtol = 1e-8

    P_oupt = 2.5e5
    f = 400.0
    TASOPT.propsys.ElectricMachine.operate_inverter!(invtr,P_oupt, f)

    @test invtr.P_input ≈ 251837.36343069648 rtol = 1e-8
end