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

    @test length(csv1[1]) == 803 # = entries w/ full ac `struct` and in default_output_indices
    @test length(csv2[1]) == 1153 # = entries w/ full ac `struct` and all output_indices

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

end #testset "io"

@testset "Composite tank input parsing" begin
    CompositeLaminate = TASOPT.materials.CompositeLaminate
    # Reads example/cryo_input.toml with optional text edits applied to the temp copy only.
    function read_cryo_variant(; extra::String = "", replacements::Vector{Pair{String,String}} = Pair{String,String}[])
        text = read(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"), String)
        for (old, new) in replacements
            @assert count(old, text) == 1
            text = replace(text, old => new)
        end
        @assert count("[Fuel.Storage]\n", text) == 1
        text = replace(text, "[Fuel.Storage]\n" => "[Fuel.Storage]\n" * extra)
        path = tempname() * ".toml"
        write(path, text)
        try
            return read_aircraft_model(path)
        finally
            rm(path, force = true)
        end
    end
    # Returns the error raised by reading a variant, or nothing if reading succeeded.
    function read_error(; kwargs...)
        try
            read_cryo_variant(; kwargs...)
            return nothing
        catch err
            return err
        end
    end

    @testset "unmodified cryo_input.toml keeps alloy walls and inert composite defaults" begin
        ac = read_cryo_variant()
        d = TASOPT.fuselage_tank()
        ft = ac.fuse_tank
        @test ft.inner_material isa StructuralAlloy
        @test ft.inner_material.name == "Al-2219-T87"
        @test ft.design_strain_limit == d.design_strain_limit
        @test ft.f_joint_add == d.f_joint_add
        @test ft.eta_composite == d.eta_composite
        @test ft.has_permeation_barrier == d.has_permeation_barrier
        @test ft.t_barrier == d.t_barrier
        @test ft.dome_layup_efficiency == d.dome_layup_efficiency
        @test ft.buckling_knockdown == d.buckling_knockdown
        @test ft.composite_min_stacks == d.composite_min_stacks
        @test ft.wall_architecture == d.wall_architecture
        @test ft.stiffener_material.name == "Al-2219-T87"
    end

    @testset "composite inner wall from explicit type and optional design keys" begin
        extra = """
            inner_vessel_material_type = "composite"
            composite_design_strain = 0.0025
            composite_joint_mass_fraction = 0.05
            composite_efficiency = 0.9
            has_permeation_barrier = true
            permeation_barrier_thickness = "0.38 mm"
            dome_layup_efficiency = 0.9
            composite_min_stacks = 2
            stiffener_material = "Ti-6242"
            """
        ac = read_cryo_variant(extra = extra,
            replacements = ["inner_vessel_material = \"Al-2219-T87\"" => "inner_vessel_material = \"CFRP-IM7-8552-PW-0090\""])
        ft = ac.fuse_tank
        @test ft.inner_material isa CompositeLaminate
        @test ft.inner_material.name == "CFRP-IM7-8552-PW-0090"
        @test ft.design_strain_limit ≈ 0.0025
        @test ft.f_joint_add ≈ 0.05
        @test ft.eta_composite ≈ 0.9
        @test ft.has_permeation_barrier == true
        @test ft.t_barrier ≈ 0.38e-3
        @test ft.dome_layup_efficiency ≈ 0.9
        @test ft.composite_min_stacks == 2
        @test ft.stiffener_material isa StructuralAlloy
        @test ft.stiffener_material.name == "Ti-6242"
    end

    @testset "composite type selection: name convenience and case-insensitive type" begin
        repl = ["inner_vessel_material = \"Al-2219-T87\"" => "inner_vessel_material = \"CFRP-IM7-8552-PW-0090\""]
        ac_name = read_cryo_variant(replacements = repl)
        @test ac_name.fuse_tank.inner_material isa CompositeLaminate
        ac_case = read_cryo_variant(extra = "inner_vessel_material_type = \"Composite\"\n", replacements = repl)
        @test ac_case.fuse_tank.inner_material isa CompositeLaminate
        @test TASOPT.build_tank_wall_material("Al-2219-T87", "ALLOY") isa StructuralAlloy
        @test TASOPT.build_tank_wall_material("CFRP-IM7-8552-PW-0090", "alloy") isa CompositeLaminate
    end

    @testset "double-walled (vacuum) tank: composite outer wall and buckling knockdown" begin
        extra = """
            outer_vessel_material_type = "composite"
            buckling_knockdown = 0.65
            """
        ac = read_cryo_variant(extra = extra, replacements = [
            "insulation_material = [\"polyurethane32\", \"polyurethane32\", \"polyurethane32\"]" =>
                "insulation_material = [\"polyurethane32\", \"vacuum\", \"polyurethane32\"]",
            "outer_vessel_material = \"Al-2219-T87\"" => "outer_vessel_material = \"CFRP-IM7-8552-PW-0090-VDA\"",
        ])
        @test ac.fuse_tank.outer_material isa CompositeLaminate
        @test ac.fuse_tank.outer_material.name == "CFRP-IM7-8552-PW-0090-VDA"
        @test ac.fuse_tank.buckling_knockdown ≈ 0.65
        @test ac.fuse_tank.inner_material isa StructuralAlloy
    end

    @testset "invalid composite inputs raise helpful errors" begin
        err = try
            TASOPT.build_tank_wall_material("Al-2219-T87", "plastic")
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("plastic", sprint(showerror, err))
        @test occursin("alloy", sprint(showerror, err))
        @test occursin("composite", sprint(showerror, err))

        @test_throws Exception TASOPT.build_tank_wall_material("Al-2219-T87", "composite")

        err_arch = read_error(extra = "wall_architecture = \"sandwich\"\n")
        @test err_arch isa ErrorException
        @test occursin("monocoque", sprint(showerror, err_arch))
    end
end #testset "Composite tank input parsing"
