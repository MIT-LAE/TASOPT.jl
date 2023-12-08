using Documenter, TASOPT

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
        "Loading and running"=>"examples/loadingrunning.md",
        "Payload-range diagram" => "examples/payload_range.md",
        "Multivariable optimization"=>"examples/optimization.md"
        ],
    "Aerodynamics" => Any[
        "aero/geometry.md",
        "aero/lift.md",
        "aero/drag.md",
        "aero/moment.md",
        "atmos/atmos.md"
        ],
    "Structures" => Any["structures/wing.md",
        "structures/fuselage.md",
        "structures/fueltanks.md"
        ],
    "Propulsion systems" => Any[
            "propulsion/propsys.md",
            "propulsion/hxfun.md",
            "propulsion/gascalc.md"
        ],
    "Stability" => "balance/balance.md",
    "Mission and sizing" => Any[
        "sizing/sizing.md",
        "sizing/weightmodels.md"
        ],
    "Miscellaneous" => Any[
        "misc/structs.md",
        "misc/dreladocs.md",
        "misc/misc.md",
        "misc/fordevs.md"
        ]
    ],
    )

deploydocs(
    repo = "github.com/MIT-LAE/TASOPT.jl.git",
)
