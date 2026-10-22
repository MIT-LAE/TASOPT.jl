# Run from the repository root: julia --project example/validation/generate_plots.jl
using TASOPT, Plots, TOML
include(TASOPT.__TASOPTindices__)

output = joinpath(@__DIR__, "generated")
mkpath(output)
for class in ("narrow", "regional", "wide")
    ac = read_aircraft_model(joinpath(@__DIR__, "../defaults/default_$class.toml"))
    @info "Sizing $class model"
    sizing = size_aircraft!(ac; iter=120, printiter=false)
    record = Dict{String,Any}("sizing_success" => sizing.success)
    try
        sizing.success || error("Aircraft sizing failed: $sizing")
        @info "Searching payload-range boundary for $class"
        pr = TASOPT.payload_range_data(ac)
        TASOPT._payload_range_plot(ac, pr; filename=joinpath(output, "$class.png"))
        merge!(record, Dict("range_nmi" => pr.range ./ 1852,
                            "payload_kg" => pr.payload ./ 9.81,
                            "initial_altitude_ft" => pr.altitude ./ 0.3048,
                            "design_feasible" => pr.design_feasible,
                            "design_scaled_residuals" => Dict(string(k)=>v for (k,v) in pairs(pr.design_report.scaled)),
                            "constraint_tolerance" => 1e-6))
    catch err
        err isa InterruptException && rethrow()
        record["error"] = sprint(showerror, err)
        @warn "No accepted diagram for $class" exception=err
    end
    @info "Writing $class results to $output"
    open(joinpath(output, "$class.toml"), "w") do io
        TOML.print(io, record; sorted=true)
    end
end
