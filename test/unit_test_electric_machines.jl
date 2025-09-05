@testset "Electric machines" begin
    #Motor
    mot = TASOPT.propsys.ElectricMachine.PMSM_motor()
    mot.N_pole_pairs = 10
    mot.J_max = 10e6
    mot.N_inverters = 10
    mot.windings.kpf = 0.35
    mot.U_max = 200.0
    mot.N_inverters = 1

    Pmech = 1e6
    shaft_speed = 10e3 #RPM
    TASOPT.propsys.ElectricMachine.size_PMSM!(mot, shaft_speed, Pmech)

    @test mot.mass ≈ 99.93627515592064 rtol = 1e-8
    @test mot.l ≈ 0.15319344713958763 rtol = 1e-8
    @test mot.P_input ≈ 1.0343052741850604e6 rtol = 1e-8
    @test mot.phase_resistance ≈ 0.0010070190689972127 rtol = 1e-8

    Pmech = 2.5e5
    shaft_speed = 5e3
    TASOPT.propsys.ElectricMachine.operate_PMSM!(mot, shaft_speed, Pmech)
    @test mot.P_input ≈ 259340.74513568674 rtol = 1e-8

    #Generator
    gen = TASOPT.propsys.ElectricMachine.PMSM_generator()
    gen.N_pole_pairs = 10
    gen.J_max = 10e6
    gen.windings.kpf = 0.35
    gen.U_max = 200.0

    Pmech = 1e6
    shaft_speed = 10e3 #RPM
    TASOPT.propsys.ElectricMachine.size_PMSM!(gen, shaft_speed, Pmech)

    @test gen.mass ≈ 99.93627515592064 rtol = 1e-8
    @test gen.l ≈ 0.15319344713958763 rtol = 1e-8
    @test gen.P_output ≈ 965694.7258149397 rtol = 1e-8
    @test gen.phase_resistance ≈ 0.0010070190689972127 rtol = 1e-8

    Pmech = 2.5e5
    shaft_speed = 5e3
    TASOPT.propsys.ElectricMachine.operate_PMSM!(gen, shaft_speed, Pmech)
    @test gen.P_output ≈ 240659.25486431326 rtol = 1e-8

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

    #Cable
    c = TASOPT.propsys.ElectricMachine.Cable()
    P = 1e5
    l = 10.0
    V = 200.0
    eff_f = c(P, V, l)

    @test eff_f(1e5) ≈ 0.9952918084 rtol = 1e-8
end