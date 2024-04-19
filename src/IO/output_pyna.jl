using NPZ, JSON, CSV, NLopt
export output_pyna, generate_pyna_dirs
"""


generates a pyNA case in its preferred location. The case comprises a 
case directory (with subdirectories as shown in `generate_pyna_dirs`)
and the necessary data to run a noise assessment.

"""
function output_pyna(pyna_path::String, ac::TASOPT.aircraft=TASOPT.load_default_model();
    case_name::String="tasopt_default", overwrite::Bool = false)

  #aircraft checks
    #ensure aircraft is sized
    if !ac.sized[1]
        throw(ArgumentError("Aircraft model provided is flagged unsized. Size and try pyna output again: "))
    end
    
    #generate the directory
    generate_pyna_dirs(pyna_path, case_name, overwrite = overwrite)

    #generate the aircraft def'n .json
    generate_pyna_json(pyna_path, ac, case_name = case_name)

    #generate the aero deck precursors and save to aircraft
    generate_pyna_aero_inputs(pyna_path, ac, case_name = case_name)

    #generate the engine deck and save to engine
    generate_pyna_engine_inputs(pyna_path, ac, case_name = case_name)

    #hope for the best

end #function


"""


generates a directory structure at the specified pyna_path*"/cases" with 
for a specific case_name. considers an overwrite flag.

pyNA/cases
└── tasopt_default
    ├── aircraft
    ├── engine
    ├── output
    ├── shielding
    └── trajectory

"""
function generate_pyna_dirs(pyna_path::String, case_name::String="tasopt_default";
    overwrite::Bool = false)

    #check if base pyna path exists
    if !isdir(pyna_path)
        throw(ArgumentError("Base pyNA path could not be found: $pyna_path"))
    end

    # Create the full path to the case directory
    case_path = joinpath(pyna_path, "cases", case_name)

    # Check if the case directory already exists
    if isdir(case_path)
        if overwrite
            # If overwrite is true, delete the existing directory
            rm(case_path; force=true, recursive=true)
        else
            # If overwrite is false, throw an error
            throw(ArgumentError("Case directory $case_name already exists, and overwrite specified false."))
        end
    end

    # Create the case directory
    mkdir(case_path)

    # Create subdirectories within the case directory
    subdirectories = ["aircraft", "engine", "output", "shielding", "trajectory"]
    for subdir in subdirectories
        mkdir(joinpath(case_path, subdir))
    end
end #function


"""

pulls aerodynamic data from sized `aircraft` model, outputs for pyna consumption 
"""
function generate_pyna_aero_inputs(pyna_path::String, ac::TASOPT.aircraft; case_name::String="tasopt_default")
    #need 3D output for pyna consistency
    #currently, only taking the default CL-CDpoint (no "alpha" variation)

    cl_array = reshape([ac.para[iaCL,iptakeoff,1]],1,1,1)
    clmax_array = cl_array #TODO: consider if this is needed at all
    cd_array = reshape([ac.para[iaCD,iprotate,1]],1,1,1)

    cl_filepath = joinpath(pyna_path, "cases", case_name, "aircraft",
                            "c_l_"*case_name*".npy")
    clmax_filepath = joinpath(pyna_path, "cases", case_name, "aircraft",
                            "c_l_max_"*case_name*".npy")
    cd_filepath = joinpath(pyna_path, "cases", case_name, "aircraft",
                            "c_d_"*case_name*".npy")
    
    npzwrite(cl_filepath, cl_array)
    npzwrite(clmax_filepath, clmax_array)
    npzwrite(cd_filepath, cd_array)
end #function


"""

pulls geometry and other data from aircraft model (able to add reasonable defaults where specific data absent)
"""
function generate_pyna_json(pyna_path::String, ac::TASOPT.aircraft; case_name::String="tasopt_default",
                            add_defaults::Bool = false)
    #initialize output dict and populate with TASOPT-sourced values, 
    #marking `nothing` if unavail (`nothing`s will be substituted with default values below)
    ac_dict = Dict( #TODO: update these. for now, test
        "mtow" => 55000.0,
        "n_eng" => nothing,
        "comp_lst" => ["wing", "tail_h", "les", "tef", "lg"],
        "af_S_h" => 20.16,
        "af_S_v" => 21.3677,
        "af_S_w" => 150.41,
        "af_b_f" => 6.096,
        "af_b_h" => 5.6388,
        "af_b_v" => 4.7244,
        "af_b_w" => 20.51304,
        "af_S_f" => 11.1484,
        "af_s" => 1.0,
        "af_d_mg" => 0.9144,
        "af_d_ng" => 0.82296,
        "af_l_mg" => 2.286,
        "af_l_ng" => 1.8288,
        "af_n_mg" => 4.0,
        "af_n_ng" => 2.0,
        "af_N_mg" => 2.0,
        "af_N_ng" => 1.0,
        "c_d_g" => 0.0240,
        "mu_r" => 0.0175,
        "B_fan" => 25,
        "V_fan" => 48,
        "RSS_fan" => 300.0,
        "M_d_fan" => 1.68,
        "inc_F_n" => 0.25,
        "TS_lower" => 0.65,
        "TS_upper" => 1.0,
        "af_clean_w" => true,
        "af_clean_h" => false,
        "af_clean_v" => true,
        "af_delta_wing" => true,
        "alpha_0" => 0
    )

    #pull sample file for defaults
    default_filepath = joinpath(__TASOPTroot__,"IO","IO_samples","pyna_aircraft_defaults_stca.json")
    def_dict = JSON.parsefile(default_filepath)
    
    #apply defaults for unspecified params
    for (key, value) in pairs(ac_dict)
        if isnothing(value)
            ac_dict[key] = def_dict[key]
        end
    end

    #output to JSON
    json_filepath = joinpath(pyna_path, "cases", case_name, "aircraft", case_name*".json")
    open(json_filepath, "w") do io
        JSON.print(io, ac_dict, 1)
    end
end #function



"""

pulls engine data from sized `aircraft` model, outputs for pyna consumption 
"""
function generate_pyna_engine_inputs(pyna_path::String, ac::TASOPT.aircraft; case_name::String="tasopt_default")

    #determine range of parameters
    alt_sample = collect(0:2000:4500) #in meters
    M_sample = collect(0:0.25:0.5)   #[-]
    TS_sample = collect(0.7:0.05:1.0)#T/Tmax (can't do past 1.0)
    # Tt4_sample = collect(1500:100:2000) #in K, proxy for thrust setting (TS) bc tasopt doesn't size that way
    # Tt4_sample = [1800, 1800]


    ip = iptakeoff  #index of flight point, doesn't do anything tfcalc! but selects the mission point within para and pare
        #TODO^check w default input to see if there's a major distinction in the first bunch
    im = 1
    icall = 2       #flag to specify thrust setting spec. approach (1 for Tt4 set, Fe computed; 2 for vice-versa) 
    icool = 2       #flag to use turbine cooling flow, 2 uses T_max_metal to set cooling fraction, fc
    initeng = 0     #flag to use current eng. vars for initialization


    pari = deepcopy(ac.pari)
    parg = deepcopy(ac.parg)
    para = deepcopy(ac.para[:,ip,im])
    pare = deepcopy(ac.pare[:,ip,im])

    data_column_labels = ["z [m]",
    "M_0 [-]",	
    "T/TMAX [-]", # = Nfb (in tasopt)
    "Fn [N]",
    "Wf [kg/s]",
    "jet_V [m/s]",
    "jet_Tt [K]",
    "jet_rho [kg/m3]",
    "jet_A [m2]",
    "jet_M [-]",
    "core_mdot_in [kg/s]",
    "core_Tt_in [K]",
    "core_Tt_out [K]",
    "core_Pt_in [Pa]",
    "core_DT_t [K]",
    "core_LPT_rho_out [kg/m3]",
    "core_HPT_rho_in [kg/m3]",
    "core_LPT_c_out [m/s]",
    "core_HPT_c_in [m/s]",
    "fan_DTt [K]",
    "fan_mdot_in [kg/s]",
    "fan_N [rpm]",
    "fan_A [m2]",
    "fan_d [m]",
    "fan_M_d [-]"
    ]

    n_runs = length(alt_sample)*length(M_sample)*length(TS_sample)
    csv_out = zeros(n_runs, length(data_column_labels))

    ctr = 1
    for alt in alt_sample
        #generate input params from the sample vars
        altkm = alt/1000. #in km
        T0, P0, _,  a0 = atmos(altkm/1000)

        #propagate input params within the par arrays
        pare[iep0] = P0
        pare[ieT0] = T0
        pare[iea0] = a0

        println("============================")
        println("alt = ", alt, " p0 = ", P0, " T0 = ", T0)
        
        for M0 in M_sample
            #propagate input within the par arrays
            pare[ieM0] = M0

            println("----------------------------")
            println("M0 = ", M0)

            #find maximum Fe at these conditions (for TS normalization)

            





            #Set function to minimize, 
            #x will be Tt4, output will be Fe
            obj(x, grad) =  Fe_tcalc_eval(x, pari, parg, para, pare, ip, 1, icool, 0) #Minimize objective function
            #icall = 1, solving for Fe from Tt4^; initeng = 0 (fresh initialization every time)

            #Use NLopt.jl to minimize function
            opt = Opt(:LN_NELDERMEAD, 1) #optimizer algorithm, dimensionality of space
            opt.lower_bounds = 700
            opt.upper_bounds = 2500
            opt.ftol_rel = 1e-3
            opt.maxeval = 100  # Set the maximum number of function evaluations

            opt.min_objective = obj
            initial_x = [1900] #K, initial estimate for max thrust Tt4

            (min_nFe,Tt4_star,returncode) = NLopt.optimize(opt, initial_x)

            if string(returncode) in ["NLOPT_MAXEVAL_REACHED", "NLOPT_FAILURE","NLOPT_OUT_OF_MEMORY","NLOPT_ROUNDOFF_LIMITED","NLOPT_FORCED_STOP"] 
                throw(ErrorException("Engine calcs for pyNA export failed to find max thrust setting, return code: "*string(returncode)))
                
            else
                Fe_max = -1*min_nFe
                println("Solve for max Fe complete: Fe_max = ", Fe_max, " @ ", Tt4_star, " K")
            end
            

            for TS in TS_sample
                println("TS = ", TS)

                pare[ieFe] = Fe_max*TS
                #does this need to be set? check in the output
                # u0 = pare[ieu0] #checked, no. leave for now as an example
            
                #execute calculation / pull the lever kronk
                tfcalc!(pari, parg, para, pare, ip, icall, icool, initeng)

                #back out what TS would be
                Tt2 = pare[ieTt2]           #fan face stag temp
                Nf = pare[ieNf]             #fan speed, dimensional
                N1 = Nf * parg[igGearf]     #LPC speed, dimensional
                N1c2 = N1 / sqrt(Tt2/Tref)  #LPC speed, corrected + dimensional

                N1des = pare[ieNbfD]        #design condition, fan speed, corrected + dimensional 
                # TS = N1c2/N1des #THIS IS NOT RIGHT. TEMP. TODO: find best way to ID the max thrust condition and speed
                
                #calculate relevant outputs
                rho8 = pare[iep8]/pare[ieR8]/pare[ieT8]
                gam8 = 1/(1-pare[ieR8]/pare[iecp8])
                M8 = pare[ieu8]/sqrt(gam8*pare[ieR8]*pare[ieT8])
                DTt_turb = pare[ieTt41] - pare[ieTt5]
                DTt_fan = pare[ieTt21] - pare[ieTt0]
                mdot_fan = pare[iemcore]*(1+pare[ieBPR])
                
                #LPT exit quantities: total density, gamma, "total speed of sound"
                rhot_49 = pare[iept49]/pare[ieTt49]/pare[ieRt49]
                gam_49 = 1/(1-pare[ieRt49]/pare[iecpt49])
                at_49 = sqrt(gam_49*pare[iecpt49]*pare[ieTt49])

                #HPT inlet quantities: ( " )
                rhot_41 = pare[iept41]/pare[ieTt41]/pare[ieRt41]
                gam_41 = 1/(1-pare[ieRt41]/pare[iecpt41])
                at_41 = sqrt(gam_41*pare[iecpt41]*pare[ieTt41])
                
                #compile outputs into the csv_row
                csv_row = [
                    alt,                #"z [m]",
                    M0,                #"M_0 [-]",	
                    TS,                #"T/TMAX [-]", # = Nfb (in tasopt)
                    pare[ieFe],             #"Fn [N]", total thrust
                    pare[ieff]*pare[iemcore],               #"Wf [kg/s]", fuel flow
                    pare[ieu8],               #"jet_V [m/s]",
                    pare[ieTt7],               #"jet_Tt [K]",
                    rho8,           #"jet_rho [kg/m3]",
                    pare[ieA8],               #"jet_A [m2]",
                    M8,               #"jet_M [-]",
                    pare[iemcore] ,               #"core_mdot_in [kg/s]",
                    pare[ieTt21] ,               #"core_Tt_in [K]",
                    pare[ieTt5] ,               #"core_Tt_out [K]",
                    pare[iept21],               #"core_Pt_in [Pa]",
                    DTt_turb ,               #"core_DT_t [K]",
                    rhot_49 ,               #"core_LPT_rho_out [kg/m3]",
                    rhot_41 ,               #"core_HPT_rho_in [kg/m3]",
                    at_49 ,               #"core_LPT_c_out [m/s]",
                    at_41 ,               #"core_HPT_c_in [m/s]",
                    DTt_fan ,               #"fan_DTt [K]",
                    mdot_fan ,               #"fan_mdot_in [kg/s]",
                    pare[ieNf]*60 ,               #"fan_N [rpm]",
                    pare[ieA2] ,               #"fan_A [m2]",
                    sqrt(pare[ieA2]/pi)*2,               #"fan_d [m]",
                    0                 #"fan_M_d [-]", 0 for now, unclear where used
                ]

                csv_out[ctr,:] = csv_row
                ctr += 1
            
            end
        end
    end

    csv_filepath = joinpath(pyna_path, "cases", case_name, "engine", "engine_deck_"*case_name*".csv")

    open(csv_filepath, "w") do io
        CSV.write(io, Tables.table(csv_out), header = data_column_labels)
    end
        #Tables.table(csv_out)
end


"""

helper function for engine input NLOpt execution
"""
function Fe_tcalc_eval(x, pari, parg, para, pare, ip, icall, icool, initeng)
    pare[ieTt4] = x[1]
    tfcalc!(pari, parg, para, pare, ip, icall, icool, initeng)
    
    nFe = -1*pare[ieFe]
    return nFe #objective fxn
end
