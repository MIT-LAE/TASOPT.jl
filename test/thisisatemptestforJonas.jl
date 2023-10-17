#Test for save_model() functionality

using TASOPT
# include("./src/TASOPT.jl")

#Part A: run the baseline default aircraft ==============
#load default
ac = load_default_model()
#run and save total weight
size_aircraft!(ac, iter=35, initwgt=false, Ldebug=false, printiter=false, saveOD=false)
igWMTO = 4
println(["Baseline MTOM: ",ac.parg[igWMTO]])

#Part B: generate, save, load, and run the tweaked aircraft ==============
#reload default and tweak max span
ac2 = load_default_model()
igbmax = 247    #max span
ac2.parg[igbmax] = 45.72 #m, = 150ft
#save to new path
newfilepath = joinpath(TASOPT.__TASOPTroot__, "IO/test_output.toml")
save_model(ac2,newfilepath)

#read tweaked model at newfilepath
ac2loaded = read_aircraft_model(newfilepath)

#solve and compare total weight
println(ac2loaded)
size_aircraft!(ac2loaded, iter=35, initwgt=false, Ldebug=false, printiter=false, saveOD=false)
println(["Tweaked MTOM: ",ac2loaded.parg[igWMTO]])