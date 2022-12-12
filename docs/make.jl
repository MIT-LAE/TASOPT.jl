using Pkg
Pkg.activate("../")

push!(LOAD_PATH, "../src")

include("../src/sizing/wsize.jl")
include("../src/aero/cdsum.jl")
include("../src/aero/fusebl.jl")
include("../src/aero/blax.jl")
include("../src/aero/blsys.jl")
include("../src/aero/axisol.jl")
include("../src/aero/trefftz.jl")
include("../src/aero/wingsc.jl")
include("../src/aero/wingpo.jl")
include("../src/aero/surfcm.jl")
include("../src/aero/surfcd.jl")
include("../src/aero/airfun2.jl")

using Documenter

makedocs(
    sitename = "TAESOPT.jl documentation",
    pages = [ "Home" => "index.md", 
    "Aerodynamics" => Any[
        "aero/geometry.md",
        "aero/lift.md",
        "aero/drag.md",
        "aero/moment.md",
    ],])
