#Test for save_model() functionality
using Revise
using TASOPT
println("----------------------------------------------")

test_choice = 1

if test_choice == 1
# =======================================================
# tests for a "careful" save into a readable toml

    # Test 0: run the baseline default aircraft ==============
    #load default
    ac1 = load_default_model()
    #run and save total weight
    igWMTO = 4


    # Test A: write unsized aircraft, then read and run 
    filepath_rewrite = joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/test_output.toml")
    save_model(ac1, filepath_rewrite)

    size_aircraft!(ac1, iter=35, initwgt=false, Ldebug=false, printiter=true, saveOD=false)
    println("Baseline MTOM: ",round(ac1.parg[igWMTO]))

    ac1_reread = read_aircraft_model(filepath_rewrite)
    size_aircraft!(ac1_reread, iter=35, initwgt=false, Ldebug=false, printiter=true, saveOD=false)
    println("Re-write and -load MTOM: ", round(ac1_reread.parg[igWMTO]))

    # Test B: write sized aircraft then read (no run)
    #TODO: waiting for output writing to be done

    # Test C: tweak baseline aircraft, save, load, and run ==============
    #reload default and tweak max span
    ac2 = load_default_model()
    #TODO: look into this
    igbmax = 247    #max span
    # ac2.parg[igbmax] = 45.72 #m, = 150ft
    #save to new path
    filepath_tweaked = joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/test_output_tweaked.toml")
    save_model(ac2,filepath_tweaked)

    #read tweaked model at newfilepath
    ac2loaded = read_aircraft_model(filepath_tweaked)

    #solve and compare total weight
    println(ac2loaded)
    size_aircraft!(ac2loaded, iter=35, initwgt=false, Ldebug=false, printiter=true, saveOD=false)
    println("Tweaked MTOM: ",round(ac2loaded.parg[igWMTO]))

elseif test_choice == 2
# ========================================================
# test for saving the struct simply and dumpily

    #indices for convenience
    igWMTO = 4      #MTOM
    imWpay = 3      #Wpayload

    #Test A - quicksave default aircraft, load it, run it
    quicksave_aircraft()
    ac1 = quickload_aircraft()
    size_aircraft!(ac1, iter=35, initwgt=false, Ldebug=false, printiter=false, saveOD=false)
    # print(ac1)
    println("Baseline MTOM, quicksave roundtrip: ",round(ac1.parg[igWMTO]))
    println("Baseline payload, quicksave, roundtrip: ", round(ac1.parm[imWpay]))

    #Test B - tweak default aircraft, quicksave it, load it, run it
    ac2 = quickload_aircraft()
    ac2.parm[imWpay] = 123456 #Newtons

    #save to new path
    filepath_tweaked = joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/testoutput_quicksave.toml")
    quicksave_aircraft(ac2,filepath_tweaked)
    ac2_reloaded = quickload_aircraft(filepath_tweaked)
    
    #comparing mtoms of the tweaked aircraft and same aircraft quicksaved
    size_aircraft!(ac2, iter=35, initwgt=false, Ldebug=false, printiter=false, saveOD=false)
    println("Tweaked MTOM, no quicksave: ",round(ac2.parg[igWMTO]))
    println("Tweaked payload, no quicksave: ", round(ac2.parm[imWpay]))
    size_aircraft!(ac2_reloaded, iter=35, initwgt=false, Ldebug=false, printiter=false, saveOD=false)
    println("Tweaked MTOM, quicksave roundtrip: ",round(ac2_reloaded.parg[igWMTO]))
    println("Tweaked payload, quicksave roundtrip: ", round(ac2_reloaded.parm[imWpay]))
end