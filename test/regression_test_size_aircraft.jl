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
                @test val1 ≈ val2 
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

    test_ac_off_design(ac, 1.0800347958058685, 152152.15193363777, 763293.7221007551)
    
    @test ac.parm[imPFEI] ≈ 0.9402313710034372 rtol=1e-4
end

@testset "Wide sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_wide.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 3.0988
    

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.1702982817495666 rtol=1e-4

end

@testset "Regional sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 1.5113

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 0.8435990587821769 rtol=1e-4

end

@testset "Hydrogen sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 2.54

    size_aircraft!(ac, iter=50; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.0042440592543385 rtol=1e-4

end

@testset "Composite example input" verbose=true begin
    include(__TASOPTindices__)
    CompositeLaminate = TASOPT.materials.CompositeLaminate
    alu_path = joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml")
    cfrp_path = joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input_composite.toml")

    @testset "aluminium example is unchanged" begin
        ac = read_aircraft_model(alu_path)
        @test ac.fuse_tank.inner_material isa StructuralAlloy
        @test ac.fuse_tank.outer_material isa StructuralAlloy
        @test ac.fuse_tank.design_strain_limit == 0.0
    end

    @testset "composite example parses" begin
        ac = read_aircraft_model(cfrp_path)
        ft = ac.fuse_tank
        @test ft.fueltype == "LH2"
        @test ft.inner_material isa CompositeLaminate
        @test ft.inner_material.name == "CFRP-IM7-8552-PW-0090"
        @test ft.outer_material isa CompositeLaminate
        @test endswith(ft.outer_material.name, "-VDA")
        @test ft.has_permeation_barrier == false
        @test ft.t_barrier ≈ 0.0
        @test ft.buckling_knockdown ≈ 0.65
        @test ft.dome_layup_efficiency ≈ 0.9
        @test ft.design_strain_limit == 0.0
    end

    @testset "composite example sizes end-to-end" begin
        ac = read_aircraft_model(cfrp_path)
        size_aircraft!(ac, iter=50; printiter=false)
        ft = ac.fuse_tank
        @test isfinite(ac.parg[igWMTO]) && ac.parg[igWMTO] > 0
        @test ft.t_wall_inner > 0
        @test ft.t_wall_inner / ft.inner_material.t_ply ≈ round(ft.t_wall_inner / ft.inner_material.t_ply) atol=1e-6
        @test ft.t_wall_outer / ft.outer_material.t_ply ≈ round(ft.t_wall_outer / ft.outer_material.t_ply) atol=1e-6

        # ASSUMPTION: composite walls should not make the aircraft heavier; wide band, check the physics sign not the number
        # like-for-like reference: the same double-walled layout with aluminium walls and no barrier
        twin = read(cfrp_path, String)
        for (a, b) in ("inner_vessel_material_type = \"composite\"" => "inner_vessel_material_type = \"alloy\"",
                       "inner_vessel_material = \"CFRP-IM7-8552-PW-0090\"" => "inner_vessel_material = \"Al-2219-T87\"",
                       "outer_vessel_material_type = \"composite\"" => "outer_vessel_material_type = \"alloy\"",
                       "outer_vessel_material = \"CFRP-IM7-8552-PW-0090-VDA\"" => "outer_vessel_material = \"Al-2219-T87\"",
                       "has_permeation_barrier = false" => "has_permeation_barrier = true")
            @assert occursin(a, twin)
            twin = replace(twin, a => b)
        end
        twin_path = tempname() * ".toml"
        write(twin_path, twin)
        ac_alu = read_aircraft_model(twin_path)
        rm(twin_path, force = true)
        size_aircraft!(ac_alu, iter=50; printiter=false)
        @test 0.8 < ac.parg[igWMTO] / ac_alu.parg[igWMTO] < 1.05
    end
end #testset "Composite example input"
