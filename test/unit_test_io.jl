#tests io functionalities
@testset "io" verbose=true begin
#A: readable TOML saves
    #check that the default model is sized identically via MTOW
        #when round-tripped via model save and read
    ac_def = load_default_model()

    filepath_rewrite = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_rewrite.toml")
    save_aircraft_model(ac_def, filepath_rewrite)
    ac_reread = read_aircraft_model(filepath_rewrite)

    size_aircraft!(ac_def, Ldebug=false, printiter=false)
    size_aircraft!(ac_reread, Ldebug=false, printiter=false)

    @test ac_def.parg[igWMTO] ≈ ac_reread.parg[igWMTO]
    rm(filepath_rewrite)

    #check via MTOW that changing an important parameter survives the save
    # and changes the solution
    ac_lopay = load_default_model()
    ac_lopay.parm[imWpay] = ac_def.parm[imWperpax, 1] * 30 #thirty pax in N
    filepath_nopay = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_nopay.toml")
    save_aircraft_model(ac_lopay, filepath_nopay)

    ac_lopay_reread = read_aircraft_model(filepath_nopay)
    size_aircraft!(ac_lopay, Ldebug=false, printiter=false)
    size_aircraft!(ac_lopay_reread, Ldebug=false, printiter=false)
    
    @test ac_lopay.parg[igWMTO] ≈ ac_lopay_reread.parg[igWMTO]
    @test !(ac_lopay.parg[igWMTO] ≈ ac_def.parg[igWMTO])
    rm(filepath_nopay)

#B: quicksaves and loads
    #test that quicksave/load roundtrip default aircraft sizes identically to default load
    filepath_quick = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_quick.toml")
    quicksave_aircraft(load_default_model(), filepath_quick)
    ac_quick = quickload_aircraft(filepath_quick)
    size_aircraft!(ac_quick, Ldebug=false, printiter=false)

    @test ac_quick.parg[igWMTO] ≈ ac_def.parg[igWMTO]
    rm(filepath_quick)

    #check via MTOW that changing an important parameter survives the quicksave
    # and changes the solution
    filepath_quick_nopay = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_quick_nopay.toml")
    ac_quick.parm[imWpay] = ac_def.parm[imWperpax, 1] * 30 #thirty pax in N
    quicksave_aircraft(ac_quick, filepath_quick_nopay)

    ac_quick_nopay_reread = quickload_aircraft(filepath_quick_nopay)
    size_aircraft!(ac_quick_nopay_reread, Ldebug=false, printiter=false)
    @test ac_quick_nopay_reread.parg[igWMTO] ≈ ac_lopay.parg[igWMTO]
    rm(filepath_quick_nopay)

    #check that a quicksave written before the `infeasibilities` field existed still loads as an aircraft
    # the old layout is mimicked by an `aircraft` struct defined without that field
    ac_old = load_default_model()
    old_fields = filter(!=(:infeasibilities), collect(fieldnames(typeof(ac_old))))
    OldLayout = Module(:OldLayout)
    Core.eval(OldLayout, Expr(:struct, true, :aircraft,
        Expr(:block, [Expr(:(::), f, fieldtype(typeof(ac_old), f)) for f in old_fields]...)))
    ac_old_layout = Core.eval(OldLayout, Expr(:call, :aircraft, [getfield(ac_old, f) for f in old_fields]...))
    filepath_quick_old = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_quick_old.jld2")
    TASOPT.JLD2.jldopen(filepath_quick_old, "w") do file
        file["ac"] = ac_old_layout
    end

    ac_quick_old = quickload_aircraft(filepath_quick_old)
    @test ac_quick_old isa TASOPT.aircraft
    @test ac_quick_old.infeasibilities == Symbol[]
    @test ac_quick_old.parm[imWpay] == ac_old.parm[imWpay]
    @test ac_quick_old.wing.layout.S == ac_old.wing.layout.S
    rm(filepath_quick_old)

#C: outputs to .csv
    using CSV
    #generate file paths
    filepath_csv = joinpath(TASOPT.__TASOPTroot__, "../test/iotest_def.csv")
    suffix = "_1"
    insertpos = length(filepath_csv) - 3
    filepath_csv2 = filepath_csv[1:insertpos-1] * suffix * filepath_csv[insertpos:end]

    #cleanup in case test was interrupted before
    isfile(filepath_csv) ? rm(filepath_csv) : nothing
    isfile(filepath_csv2) ? rm(filepath_csv2) : nothing
    
    #generate output files, indices are default_output_indices
    output_csv(ac_def, filepath_csv)
    output_csv(ac_def, filepath_csv, 
                includeMissions = true, includeFlightPoints = true)
    output_csv(ac_def, filepath_csv,
                includeMissions = false, includeFlightPoints = true,
                forceMatrices = true)
    output_csv(ac_def, filepath_csv, includeFlightPoints = true)

    #this call generates the second file since indices don't match and fuse_tank excluded
    output_csv(ac_def, filepath_csv, indices = Dict(), struct_excludes = ["fuse_tank"])
    
    #test that it creates the files as expected
    @test isfile(filepath_csv)
    @test isfile(filepath_csv2)
    
    #pull the files
    csv1 = CSV.File(filepath_csv)
    csv2 = CSV.File(filepath_csv2)

    #check row and column counts
    @test size(csv1,1) == 4 #4 rows w default indices
    @test size(csv2,1) == 1 #1 row with addl indices (all)

    @test length(csv1[1]) == 987 # = entries w/ full ac `struct` and in default_output_indices
    @test length(csv2[1]) == 1351 # = entries w/ full ac `struct` and all output_indices

    #test the nested vectors within par arrays
    #a: row 1 in both csvs matches the design cruise point/mission 
    @test parse(Float64, string(csv1[1].iaalt)) == ac_def.para[iaalt,ipcruise1,1]
    @test parse(Float64, string(csv2[1].iaalt)) == ac_def.para[iaalt,ipcruise1,1]
    
    #b: row 2 has the correct structure (m flight points, n missions)
    #note - for simplicity of imports, evaluate structure by counting brackets
    # since much of the data is parsed as Strings when using CSV.File
    #if more than one mission, one more bracket added
    entry2 = csv1[2].iaalt
    bracket_for_missions = Int( (size(ac_def.parm,2) > 1) ) 
    @test count(ichar->(ichar=='['), entry2) == iptotal + bracket_for_missions
    #c: row 3 has the same structure as b despite 1 mission bc forceMatrices
    entry3 = csv1[3].iaalt
    @test count(ichar->(ichar=='['), entry3) == iptotal + bracket_for_missions
    #d: row 4 only has flight points, thus 1 bracket
    entry4 = csv1[4].iaalt
    @test count(ichar->(ichar=='['), entry4) == 1

    #test that structs are properly saved
    #a: check fuselage detail, wing detail
    @test csv1[1][Symbol("fuselage.layout.cross_section.radius")] == ac_def.fuselage.layout.cross_section.radius
    @test csv2[1][Symbol("wing.inboard.cross_section.thickness_to_chord")] == ac_def.wing.inboard.cross_section.thickness_to_chord

    #b: test exclusion of structs when indicated
    @test !haskey(csv2[1], Symbol("fuse_tank.tank_type")) #excluded

    #cleanup
    rm(filepath_csv)
    rm(filepath_csv2)

#D: wing strut material is read from the input file
    @testset "strut material" begin
        al_default = StructuralAlloy("TASOPT-Al")
        same_alloy(a, b) = all(getfield(a, f) == getfield(b, f) for f in (:name, :ρ, :E, :σmax, :τmax))
        read_toml_edit(edit) = mktempdir() do dir
            path = joinpath(dir, "strut_material_test.toml")
            write(path, edit(read(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml"), String)))
            read_aircraft_model(path)
        end

        #default keys leave the strut material unchanged
        @test same_alloy(ac_def.wing.strut.caps.material, al_default)

        #changed keys build the strut caps material, without touching the wing caps
        ac_strut = read_toml_edit(txt -> replace(txt,
            "strut_material = \"TASOPT-Al\"" => "strut_material = \"Al-7075\"",
            "strut_max_avg_stress = 1.1" => "strut_max_avg_stress = 1.2",
            "strut_safety_factor = 1.5" => "strut_safety_factor = 2.0"))
        al7075 = StructuralAlloy("Al-7075", max_avg_stress = 1.2, safety_factor = 2.0)
        @test same_alloy(ac_strut.wing.strut.caps.material, al7075)
        @test ac_strut.wing.strut.caps.material.σmax ≈ al7075.YTS / 1.2 / 2.0
        @test same_alloy(ac_strut.wing.inboard.caps.material, ac_def.wing.inboard.caps.material)

        #legacy strut stress key is ignored
        ac_legacy = read_toml_edit(txt -> replace(txt, "sigma_struts  = \"30e3 lbf/in2\"" => "sigma_struts  = \"10e3 lbf/in2\""))
        @test ac_legacy.wing.strut.caps.material.σmax == ac_def.wing.strut.caps.material.σmax
    end

end #testset "io"
