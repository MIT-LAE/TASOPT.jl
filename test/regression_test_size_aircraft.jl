# Define a function to check if each value in two structs is equal
function check_struct_equivalence(s1, s2)
    fields_s1 = fieldnames(typeof(s1))
    fields_s2 = fieldnames(typeof(s2))
    
    # Check if both structs have the same fields
    if fields_s1 != fields_s2
        return false
    end
    
    # Check if each field has the same value in both structs
    for field in fields_s1
        val1 = getproperty(s1, field)
        val2 = getproperty(s2, field)
        if typeof(val1) == typeof(val2)
            if typeof(val1) != Float64
                if !check_struct_equivalence(val1, val2)
                    return false
                end
            else
                # println(field)
                @test (isnan(val1) && isnan(val2)) || val1 ≈ val2 
            end
        else
            return false
        end
    end
    
    return true
end

#Simple function to call fly_mission!() and test on- and off-design performance
function test_ac_off_design(ac, PFEI, Wfuel, WTO)
    @testset "Off-design" begin
        TASOPT.fly_mission!(ac, 2; printTO=false)

        @test ac.parm[imPFEI, 2] ≈ PFEI
        @test ac.parm[imWfuel, 2] ≈ Wfuel
        @test ac.parm[imWTO, 2] ≈ WTO
    end
end

@testset "Default sizing" verbose=true begin
    ac = load_default_model()
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius == 1.9558
    
    include(joinpath(TASOPT.__TASOPTroot__, "../test/default_sized.jl"))
    # Fuselage
    include(joinpath(TASOPT.__TASOPTroot__, "../test/default_structures.jl"))

    size_aircraft!(ac; printiter=false);

    @testset "Fuselage" begin
        @test  check_struct_equivalence(ac_test.fuselage, ac.fuselage)
    end

    @testset "Wing" begin
        @test  check_struct_equivalence(ac_test.wing, ac.wing)
    end

    @testset "Htail" begin
        @test  check_struct_equivalence(ac_test.htail, ac.htail)
    end

    @testset "Vtail" begin
        @test  check_struct_equivalence(ac_test.vtail, ac.vtail)
    end

    @testset "Geometry" begin
        for i in eachindex(parg)
            @test parg[i] ≈ ac.parg[i] 
        end
    end

    @testset "Aero" begin
        for i in eachindex(para)
            @test para[i] ≈ ac.para[i] 
        end
    end

    @testset "Propulsion" begin
        for i in eachindex(pare)
            @test pare[i] ≈ ac.pare[i] 
        end
    end

    test_ac_off_design(ac, 1.0800347958058674, 152152.15193363762, 763293.7221007552)
    
    @test ac.parm[imPFEI] ≈ 0.9402313710034373 rtol=1e-4
end

@testset "Wide sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_wide.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 3.0988
    

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.1702982817493792 rtol=1e-4

end

@testset "Regional sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 1.5113

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 0.843599058782176 rtol=1e-4

end

@testset "Hydrogen sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 2.54

    size_aircraft!(ac, iter=50; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.0042440508631023 rtol=1e-4

end

# Evaluates reference files from test/ in a module of their own, so that the variables they define
# (parg, parm, para, pare, ac_test, fuse, wing, ...) do not leak into the other testsets.
function include_references(files...)
    reference = Module(:Reference)
    Core.eval(reference, :(using TASOPT))
    for file in (TASOPT.__TASOPTindices__, (joinpath(TASOPT.__TASOPTroot__, "../test", f) for f in files)...)
        Base.include(reference, file)
    end
    return reference
end

# Compares the sized aircraft `ac` with the structures and parameter arrays stored in `reference`
function test_against_reference(ac, reference)
    @testset "Fuselage" begin
        @test check_struct_equivalence(reference.ac_test.fuselage, ac.fuselage)
    end

    @testset "Wing" begin
        @test check_struct_equivalence(reference.ac_test.wing, ac.wing)
    end

    @testset "Htail" begin
        @test check_struct_equivalence(reference.ac_test.htail, ac.htail)
    end

    @testset "Vtail" begin
        @test check_struct_equivalence(reference.ac_test.vtail, ac.vtail)
    end

    @testset "Geometry" begin
        for i in eachindex(reference.parg)
            @test reference.parg[i] ≈ ac.parg[i]
        end
    end

    @testset "Mission" begin
        for i in eachindex(reference.parm)
            @test reference.parm[i] ≈ ac.parm[i]
        end
    end

    @testset "Aero" begin
        for i in eachindex(reference.para)
            @test reference.para[i] ≈ ac.para[i]
        end
    end

    @testset "Propulsion" begin
        for i in eachindex(reference.pare)
            @test reference.pare[i] ≈ ac.pare[i]
        end
    end
end

# Strut sizing criterion, governing wing load cases, strut attachment reinforcement and infeasibility
# flags of a sized truss-braced aircraft. These are symbols that the reference files do not record.
function test_strut_cases(ac; strut_criterion, governing_root, governing_break,
                          governing_strut, governing_strut_compression, infeasibilities)
    @test ac.options.opt_wing_config == WingConfig.TrussBraced
    @test ac.wing.strut.sizing_criterion == strut_criterion
    @test ac.wing.sizing_loads.governing_root == governing_root
    @test ac.wing.sizing_loads.governing_break == governing_break
    @test ac.wing.sizing_loads.governing_strut == governing_strut
    @test ac.wing.sizing_loads.governing_strut_compression == governing_strut_compression
    @test ac.fuselage.strut_attach_criterion == :strut_weight
    @test ac.wing.jury.weight > 0.0 # the truss-braced wing has a jury strut
    @test Set(ac.infeasibilities) == infeasibilities
end

@testset "TBW sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_tbw.toml"))

    include(__TASOPTindices__)

    reference = include_references("default_tbw_sized.jl", "default_tbw_structures.jl")

    size_aircraft!(ac; iter=100, printiter=false)

    test_against_reference(ac, reference)

    @testset "Headline numbers" begin
        @test ac.options.opt_fuel == FuelType.JetA
        @test ac.parg[igWMTO] ≈ 809446.8407033044 rtol=1e-4
        @test ac.wing.weight ≈ 56061.07398635441 rtol=1e-4
        @test ac.wing.strut.weight ≈ 24596.220913845293 rtol=1e-4
        @test ac.parg[igWfuel] ≈ 224605.59205944135 rtol=1e-4
        @test ac.parm[imPFEI] ≈ 0.8857429625554126 rtol=1e-4
        @test ac.para[iaCL, ipcruise1, 1] / ac.para[iaCD, ipcruise1, 1] ≈ 17.430988954584222 rtol=1e-4
    end

    @testset "Strut" begin
        @test ac.wing.strut.area ≈ 0.021411926295878153 rtol=1e-4
        @test ac.wing.strut.axial_force ≈ 2.1102105409874977e6 rtol=1e-4
        @test ac.fuselage.strut_attach.W ≈ 1475.7732548307174 rtol=1e-4
        test_strut_cases(ac; strut_criterion = :buckling,
                         governing_root = :taxi, governing_break = :pullup,
                         governing_strut = :pullup, governing_strut_compression = :pushover,
                         infeasibilities = Set([:wing_fuel_volume_exceeded]))
    end

    test_ac_off_design(ac, 1.0300480637746714, 145110.1669196255, 701260.3836634884)
end

@testset "TBW LCH4 smoke" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_tbw_lch4.toml"))

    include(__TASOPTindices__)

    @test ac.options.opt_wing_config == WingConfig.TrussBraced
    @test ac.options.opt_fuel == FuelType.CH4

    size_aircraft!(ac; iter=100, printiter=false)
    @test ac.is_sized == [true]

    # Loose checks only: no golden values for this aircraft
    @test isfinite(ac.parg[igWMTO])
    @test ac.parg[igWMTO] ≈ 872e3 rtol=0.15
    for weight in (ac.wing.weight, ac.wing.strut.weight, ac.fuselage.weight, ac.parg[igWftank])
        @test isfinite(weight) && weight > 0
    end
    @test ac.wing.strut.sizing_criterion in (:tension, :buckling, :gauge)
    @test issubset(Set(ac.infeasibilities), Set([:wing_beam_column_margin, :wing_fuel_volume_exceeded]))
end
