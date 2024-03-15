#tests io functionalities

@testset "io" begin
    
#A: readable TOML saves
    @testset "toml model io" begin
    #check that the default model is sized identically via MTOW
        #when round-tripped via model save and read
    ac_def = load_default_model()

    filepath_rewrite = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_rewrite.toml")
    save_model(ac_def, filepath_rewrite)
    ac_reread = read_aircraft_model(filepath_rewrite)

    size_aircraft!(ac_def, Ldebug=false, printiter=false, saveOD=false)
    size_aircraft!(ac_reread, Ldebug=false, printiter=false, saveOD=false)

    @test ac_def.parg[igWMTO] ≈ ac_reread.parg[igWMTO]
    rm(filepath_rewrite)

    #check via MTOW that changing an important parameter survives the save
    # and changes the solution
    ac_nopay = load_default_model()
    ac_nopay.parm[imWpay] = 1 #N
    filepath_nopay = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_nopay.toml")
    save_model(ac_nopay, filepath_nopay)

    ac_nopay_reread = read_aircraft_model(filepath_nopay)
    size_aircraft!(ac_nopay, Ldebug=false, printiter=false, saveOD=false)
    size_aircraft!(ac_nopay_reread, Ldebug=false, printiter=false, saveOD=false)
    
    @test ac_nopay.parg[igWMTO] ≈ ac_nopay_reread.parg[igWMTO]
    @test !(ac_nopay.parg[igWMTO] ≈ ac_def.parg[igWMTO])
    rm(filepath_nopay)
    end #testset "toml model io"
#B: quicksaves and loads



#C: outputs to .csv





end #testset "io"




