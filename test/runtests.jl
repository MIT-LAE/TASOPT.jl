using TASOPT
using Test

@testset verbose=true "TASOPT" begin
    include("unit_test_structures.jl")
    include("unit_test_aero.jl")
    include("regression_test_wsize.jl")
end