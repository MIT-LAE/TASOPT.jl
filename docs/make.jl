using Pkg
Pkg.activate("../")

push!(LOAD_PATH, "../src")

include("../src/TASOPT.jl")
aerodynamics = TASOPT.aerodynamics
structures = TASOPT.structures
aircraft = TASOPT.aircraft

using Documenter

makedocs(
    sitename = "TASOPT.jl documentation",
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
    ],
    "Stability" => "balance/balance.md",
    "Sizing" => "sizing/sizing.md",
    "Miscellaneous" => Any[
        "misc/aircraft.md",
        "misc/misc.md"
    ]
    ])
