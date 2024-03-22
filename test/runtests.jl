using TASOPT
using Test

@testset verbose=true "TASOPT" begin
    include("unit_test_structures.jl")
    include("unit_test_aero.jl")
    include("regression_test_wsize.jl")
    include("unit_test_heat_exchanger.jl")
    include("unit_test_PEMFC.jl")
    include("unit_test_materials.jl")
    #include("unit_test_fueltank.jl")
    include("unit_test_outputs.jl")
    include("unit_test_io.jl")

end