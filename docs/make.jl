using Pkg
Pkg.activate("../")

push!(LOAD_PATH, "../src")

include("../tasopt.jl")

using Documenter

makedocs(
    sitename = "TAESOPT.jl documentation",
    pages = [ "Home" => "index.md", 
    "Atmospheric properties" => "atmos/atmos.md",
    "Aerodynamics" => Any[
        "aero/geometry.md",
        "aero/lift.md",
        "aero/drag.md",
        "aero/moment.md",
    ],
    "Structures" => Any[
        "structures/wing.md",
        "structures/fuselage.md",
        "structures/fueltanks.md"
    ],
    "Propulsion systems" => Any[
        "propulsion/propsys.md",
    ]
    ])
