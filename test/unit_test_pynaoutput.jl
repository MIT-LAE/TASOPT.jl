




@testset "generate_pyna_dirs" begin
    #throws error if base pyna path doesn't exist
    @test_throws ArgumentError generate_pyna_dirs("nonexistent/path", "tasopt_default", false)

    # Test case where directory doesn't exist and overwrite is false
    #set up the pyna-like temp directory
    tmp_dir = mktempdir()
    mkdir(joinpath(tmp_dir, "cases"))

    generate_pyna_dirs(tmp_dir, "test_case", false)
    @test all(isdir, [joinpath(tmp_dir, "cases", "test_case", subdir) for subdir in ["aircraft", "engine", "output", "shielding", "trajectory"]])
    rm(tmp_dir; force=true, recursive=true)


    # Test case where directory already exists and overwrite is true
    tmp_dir = mktempdir()
    mkdir(joinpath(tmp_dir, "cases"))
    mkdir(joinpath(tmp_dir, "cases", "existing_case"))
    generate_pyna_dirs(tmp_dir, "existing_case", true)
    @test all(isdir, [joinpath(tmp_dir, "cases", "existing_case", subdir) for subdir in ["aircraft", "engine", "output", "shielding", "trajectory"]])
    rm(tmp_dir; force=true, recursive=true)

    # Test case where directory already exists and overwrite is false
    tmp_dir = mktempdir()
    mkdir(joinpath(tmp_dir, "cases"))
    mkdir(joinpath(tmp_dir, "cases", "existing_case"))
    @test_throws ArgumentError generate_pyna_dirs(tmp_dir, "existing_case", false)
end
