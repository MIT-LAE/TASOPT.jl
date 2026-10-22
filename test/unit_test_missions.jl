@testset "mission" verbose=true begin

    ac = load_default_model()
    size_aircraft!(ac; printiter=false);

    #a. verify that fly_mission!() produces same results as size_aircraft!() call
    pfei_sizecall = ac.parm[imPFEI, 1]

    fly_mission!(ac, 1; printTO=false)
    pfei_missioncall = ac.parm[imPFEI, 1]
    @test pfei_sizecall ≈ pfei_missioncall

    @test TASOPT.takeoff!(ac; printTO=false)

    @testset "takeoff failure reporting" begin
        ac_engine_out = deepcopy(ac)
        ac_engine_out.parg[igneng] = 1.0
        ac_engine_out.parg[igmuroll] = 0.0
        @test !TASOPT.takeoff!(ac_engine_out; printTO=false)
        @test isfinite(ac_engine_out.parm[imgamVTO, 1])
        @test isfinite(ac_engine_out.parm[imgamVBF, 1])

        @test !TASOPT.takeoff!(ac; printTO=false, itermax=0)

        balance_converged = @test_logs (:warn, r"Pitch not converged") TASOPT.size_htail(
            ac,
            view(ac.para, :, ipdescentn, 1),
            view(ac.para, :, ipcruise1, 1),
            view(ac.para, :, ipcruise1, 1);
            itermax=0,
        )
        @test !balance_converged
    end
        
end
