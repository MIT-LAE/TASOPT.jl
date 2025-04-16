using Documenter, TASOPT

push!(LOAD_PATH, "../src")

# Shorthands for convenience
aerodynamics = TASOPT.aerodynamics
structures = TASOPT.structures
engine = TASOPT.engine
airc = TASOPT.aircraft      #"aircraft = ..." conflicted with the global namespace
CryoTank = TASOPT.CryoTank
propsys = TASOPT.propsys

makedocs(
    repo = Documenter.Remotes.GitHub("MIT-LAE", "TASOPT.jl"),
    sitename = "TASOPT.jl",
    
    pages = [ "Home" => "index.md", 
    "Examples" => Any[
        "Loading and running"=>"examples/loadingrunning.md",
        "Payload-range diagram" => "examples/payload_range.md",
        "Optimization" => Any[
            "Nelder Mead"=>"examples/NM_optimization.md",
            "Gradient Based Optimization"=>"examples/gradient_based_optimization.md",
        ],
        "Sensitivities"=>"examples/sensitivity.md"
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
        "structures/landing_gear.md",
        "structures/cabin_sizing.md"
        ],
    "Propulsion systems" => Any[
            "propulsion/propsys.md",
            "propulsion/gascalc.md",
            "propulsion/hxfun.md",
            "propulsion/PEMfuelcell.md",
            "propulsion/ElectricMachines.md",
        ],
    "Cryogenic tanks" => Any[
            "cryo_tank/fueltanks.md",
            "cryo_tank/cryotank.md"
        ],
    "Stability" => "balance/balance.md",
    "Mission and sizing" => Any[
        "sizing/sizing.md",
        "sizing/weightmodels.md"
        ],
    "Data and I/O" => Any[
        "data_io/data_basics.md",
        "data_io/structs.md",
        "data_io/data_io.md"
        ],

    "Miscellaneous" => Any[
        "misc/dreladocs.md",
        # "misc/misc.md",
        "misc/fordevs.md"
        ]
    ],

    format = Documenter.HTML(; mathengine=
        Documenter.KaTeX(
            Dict(:delimiters => [
                Dict(:left => raw"$",   :right => raw"$",   display => false),
                Dict(:left => raw"$$",  :right => raw"$$",  display => true),
                Dict(:left => raw"\[",  :right => raw"\]",  display => true),
                ],
                :macros => 
                Dict("\\genfuel" =>
                        raw"{\mathrm{C}_{x_{\mathrm{C}}}\mathrm{H}_{x_{\mathrm{H}}}\mathrm{O}_{x_{\mathrm{O}}}\mathrm{N}_{x_{\mathrm{N}}}}",
                    raw"\Ru" => raw"R_{\mathrm{univ.}}",
                    raw"\Pstd" => raw"P_{\mathrm{std}}",
                    raw"\Tstd" => raw"T_{\mathrm{std}}",
                ),
            )
        )
    )

)

deploydocs(
    repo = "github.com/MIT-LAE/TASOPT.jl.git",
)
