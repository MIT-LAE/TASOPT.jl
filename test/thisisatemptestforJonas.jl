#Test for save_model() functionality

using TASOPT
# include("./src/TASOPT.jl")

# Test 0: run the baseline default aircraft ==============
#load default
ac = load_default_model()
#run and save total weight
size_aircraft!(ac, iter=35, initwgt=false, Ldebug=false, printiter=false, saveOD=false)
igWMTO = 4
println("Baseline MTOM: ",round(ac.parg[igWMTO]))

# Test A: write unsized then read and run 
filepath_rewrite = joinpath(TASOPT.__TASOPTroot__, "IO/test_output.toml")
save_model(ac, filepath_rewrite)
ac_reread = read_aircraft_model(filepath_rewrite)
size_aircraft!(ac_reread, iter=35, initwgt=false, Ldebug=false, printiter=false, saveOD=false)
println("Re-write and -load MTOM: ", round(ac.parg[igWMTO]))

# Test B: write sized aircraft then read (no run)
#TODO: waiting for output writing to be done

# Test C: tweak baseline aircraft, save, load, and run ==============
#reload default and tweak max span
ac2 = load_default_model()
igbmax = 247    #max span
ac2.parg[igbmax] = 45.72 #m, = 150ft
#save to new path
filepath_tweaked = joinpath(TASOPT.__TASOPTroot__, "IO/test_output_tweaked.toml")
save_model(ac2,filepath_tweaked)

#read tweaked model at newfilepath
ac2loaded = read_aircraft_model(filepath_tweaked)

#solve and compare total weight
println(ac2loaded)
size_aircraft!(ac2loaded, iter=35, initwgt=false, Ldebug=false, printiter=false, saveOD=false)
println("Tweaked MTOM: ",round(ac2loaded.parg[igWMTO]))