using Pkg
Pkg.activate("..")
using TASOPT

push!(LOAD_PATH, "../src")

# Shorthands for convenience
aerodynamics = TASOPT.aerodynamics
structures = TASOPT.structures
engine = TASOPT.engine
aircraft = TASOPT.aircraft

using Documenter
makedocs(
    remotes = nothing,
    sitename = "TASOPT.jl Documentation",
    pages = [ "Home" => "index.md", 
    "Examples" => Any[
        "Payload Range Diagram" => "examples/payload_range.md",
        "Multivariable Optimization"=>"examples/optimization.md"
        ],
    "Atmospheric properties" => "atmos/atmos.md",
    "Aerodynamics" => Any[
        "aero/geometry.md",
        "aero/lift.md",
        "aero/drag.md",
        "aero/moment.md",
        ],
        "Structures" => Any["structures/wing.md",
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
            "misc/misc.md",
            "misc/dreladocs.md"
        ]
    ],
    )
