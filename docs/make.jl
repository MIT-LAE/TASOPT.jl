push!(LOAD_PATH, "../src")

include("../src/sizing/wsize.jl")

using Documenter

makedocs(sitename="TAESOPT.jl document")