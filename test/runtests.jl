using TASOPT
using Test

@testset verbose=true "TASOPT" begin
    include("unit_test_structures.jl")
    include("unit_test_aero.jl")
    include("unit_test_heat_exchanger.jl")
    include("regression_test_wsize.jl")
    include("unit_test_fueltank.jl")
    include("unit_test_outputs.jl")

end