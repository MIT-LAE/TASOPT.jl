using BenchmarkTools
using Profile
using ProfileView
using TASOPT

#----Benchmark all files----
println("---------------------------------------")
println("Starting Benchmarks")
println("---------------------------------------")
println("Benchmarking DEFAULT aircraft")
ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml")) # MODIFY <path> appropriately
res1 = @benchmark size_aircraft!(ac, iter=50; printiter = false) seconds=30 samples=5

println("Benchmarking REGIONAL aircraft")
ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml")) # MODIFY <path> appropriately
res2 = @benchmark size_aircraft!(ac, iter=50; printiter = false) seconds=30 samples=5

println("Benchmarking WIDEBODY aircraft")
ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_wide.toml")) # MODIFY <path> appropriately
res3 = @benchmark size_aircraft!(ac, iter=50; printiter = false) seconds=30 samples=5

println("Benchmarking HYDROGEN aircraft")
ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml")) # MODIFY <path> appropriately
res4 = @benchmark size_aircraft!(ac, iter=50; printiter = false) seconds=30 samples=5

println("---------------------------------------")
println("Benchmark Results")
println("---------------------------------------")
println("DEFAULT aircraft")
display(res1)
println("REGIONAL aircraft")
display(res2)
println("WIDEBODY aircraft")
display(res3)
println("HYDROGEN aircraft")
display(res4)

#----Profile the code----
#Select aircraft file to run profile with
ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml")) 
Profile.clear()
@profile size_aircraft!(ac, iter=50; printiter = false)
#Profile.print() #Print the profile results (prints large set of results)

ProfileView.view() #View the profile