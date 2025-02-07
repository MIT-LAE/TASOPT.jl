""" 
    tfwrap!(ac, case, imission, ip, initialize_engine, iterw = 0)

Calls the turbofan sizing or off-design performance functions for the aircraft's
turbofan model. This function is basically a wrapper on tfcalc!, going from the
basic engine inputs to those required by the function and storing the outputs.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object
    - `case::String`: case identifier, e.g. "sizing" or "off_design"
    - `imission::Int64`: mission index
    - `ip::Int64`: mission point index
    - `initialize_engine::Bool`: flag to initialize engine 
        0  initialize variables for iteration in engine
        1  use current variables as initial guesses in engine
    - `iterw::Int64`: sizing loop iteration

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine parameters.
"""
function tfwrap!(ac, case::String, imission::Int64, ip::Int64, initialize_engine::Bool, iterw::Int64 = 0)
    #Unpack data storage arrays
    pari = ac.pari
    parg = ac.parg
    para = view(ac.para, :, :, imission)
    pare = view(ac.pare, :, :, imission)
    wing = ac.wing
    
    if case == "design"
        icall = 0
        icool = 1
        if (iterw == 1 || (~initialize_engine))
            # initialize engine state
            initialize_engine_firstiter  = 0
        else
            # start with current engine state
            initialize_engine_firstiter  = 1
        end

        ichoke5, ichoke7 = tfcalc!(pari,parg,view(para, :, ip), view(pare, :, ip), wing, ip, icall, icool, initialize_engine_firstiter )

        # store engine design-point parameters for all operating points
        parg[igA5] = pare[ieA5, ip] / pare[ieA5fac, ip]
        parg[igA7] = pare[ieA7, ip] / pare[ieA7fac, ip]
        
        pare[ieA2, :] .= pare[ieA2, ip]
        pare[ieA25, :] .= pare[ieA25, ip]
        pare[ieA5, :] .= parg[igA5] .* pare[ieA5fac, :]
        pare[ieA7, :] .= parg[igA7] .* pare[ieA7fac, :]

        pare[ieNbfD, :] .= pare[ieNbfD, ip]
        pare[ieNblcD, :] .= pare[ieNblcD, ip]
        pare[ieNbhcD, :] .= pare[ieNbhcD, ip]
        pare[ieNbhtD, :] .= pare[ieNbhtD, ip]
        pare[ieNbltD, :] .= pare[ieNbltD, ip]

        pare[iembfD, :] .= pare[iembfD, ip]
        pare[iemblcD, :] .= pare[iemblcD, ip]
        pare[iembhcD, :] .= pare[iembhcD, ip]
        pare[iembhtD, :] .= pare[iembhtD, ip]
        pare[iembltD, :] .= pare[iembltD, ip]

        pare[iepifD, :] .= pare[iepifD, ip]
        pare[iepilcD, :] .= pare[iepilcD, ip]
        pare[iepihcD, :] .= pare[iepihcD, ip]
        pare[iepihtD, :] .= pare[iepihtD, ip]
        pare[iepiltD, :] .= pare[iepiltD, ip]
        
    elseif case == "off_design"
        icall = (ip in range(ipstatic, ipclimbn)) ? 1 : 2 #One if in range, 2 if not
        icool = 1
        ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), wing, ip, icall, icool, Int(initialize_engine))

    elseif case == "cooling_sizing"
        icall = 1
        icool = 2
        ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), wing, ip, icall, icool, Int(initialize_engine))

        # Tmetal was specified... set blade row cooling flow ratios for all points
        for jp = 1:iptotal
            for icrow = 1:ncrowx
                pare[ieepsc1+icrow-1, jp] = pare[ieepsc1+icrow-1, ip]
            end
            # also set first estimate of total cooling mass flow fraction
            pare[iefc, jp] = pare[iefc, ip]
        end
    end
end