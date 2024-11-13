using TOML
export read_aircraft_model, load_default_model

"""
    read_input(k::String, dict::AbstractDict=data, 
    default_dict::AbstractDict = default)

Reads the input from a given dictonary (typically parsed from a TOML file).
If requested input does not exist in dictionary, looks for value in default input
and stores default value into the given dictionary (primarily for later output/
saving as an aircraft model file)
"""
function read_input(k::String, dict::AbstractDict=data, 
    default_dict::AbstractDict = default)

    get!(dict, k) do 
        if k in keys(default_dict)
            println("\n")
            @info """$k not found in user specified input file. 
            Reading $k from default TASOPT input:
            \n$k = $(default_dict[k])\n\n"""
            default_dict[k]
        else
            error("Requested key/parameter is not supported. Check the default 
            input file to see all available input options. Key: "*k)
        end
    end
end

function get_default_input_file(input)
    if input == 0 || lowercase(string(input)) == "regional aircraft"
        defaultfile = joinpath(TASOPT.__TASOPTroot__, "IO/default_regional.toml")
    elseif input == 1 || lowercase(string(input)) == "narrow body aircraft"
        defaultfile = joinpath(TASOPT.__TASOPTroot__, "IO/default_input.toml")
    elseif input == 2 || lowercase(string(input)) == "wide body aircraft"
        defaultfile = joinpath(TASOPT.__TASOPTroot__, "IO/default_wide.toml")
    else
        println("\n")
        @info """Requested aircraft type not supported. Selecting Narrow Body Aircraft"""
        defaultfile = joinpath(TASOPT.__TASOPTroot__, "IO/default_input.toml")
    end
    return TOML.parsefile(defaultfile)
end
# Convenience functions to convert to SI units
Speed(x)    = convertSpeed(parse_unit(x)...)
Distance(x)      = convertDist(parse_unit(x)...)
Force(x)    = convertForce(parse_unit(x)...)
Pressure(x) = convertPressure(parse_unit(x)...)
Stress = Pressure
Density(x)  = convertDensity(parse_unit(x)...)
Area(x)     = convertArea(parse_unit(x)...)
Vol(x)      = convertVolume(parse_unit(x)...)
Angle(x)    = convertAngle(parse_unit(x)...)
Time(x)     = convertTime(parse_unit(x)...)
Temp(x)     = convertTemp(parse_unit(x)...)

"""
    read_aircraft_model(datafile; 
    defaultfile = joinpath(TASOPT.__TASOPTroot__, "IO/default_input.toml"))

Reads a specified TOML file that describes a TASOPT `aircraft` model 
with a fall back to the default `aircraft` definition 
provided in \"src/IO/default_input.toml\""

!!! note "Deviating from default"
    Extending `read_input.jl` and `save_model.jl` is recommended for models deviating appreciably 
    from the default functionality. Thorough knowledge of the model is required.

# Examples
```julia-repl
julia> read_aircraft_model("src/IO/input.toml")


┌ Info: engine_location not found in user specified input file. 
│ Reading engine_location from default TASOPT input:
│ 
│ engine_location = wing
└ 

┌ Info: pylon_weight_fraction not found in user specified input file. 
│ Reading pylon_weight_fraction from default TASOPT input:
│ 
│ pylon_weight_fraction = 0.1
└ 
Name: Example TASOPT Model;
Wpay = 210.0 kN
Des. Range  = 5.56e6 km
Cruise Mach = 0.8

```
"""
function read_aircraft_model(
    datafile=joinpath(TASOPT.__TASOPTroot__, "IO/default_input.toml"); 
    defaultfile = joinpath(TASOPT.__TASOPTroot__, "IO/default_input.toml"))

data = TOML.parsefile(datafile)
default = TOML.parsefile(defaultfile)
ac_descrip = get(data, "Aircraft Description", Dict{})
name = get(ac_descrip, "name", "Untitled Model")
description = get(ac_descrip, "description", "---")
sized = get(ac_descrip, "sized",[false])

# Check if aircraft type exists
if "aircraft_type" in keys(ac_descrip)
    default = get_default_input_file(ac_descrip["aircraft_type"])
    ac_type_fixed = true
else
    println("\n")
    @info """Aircraft Type not specified: Selecting AC type depending on range or payload"""
    ac_type_fixed = false
end

#Get number of missions to create data arrays
mis = read_input("Mission", data, default)
dmis = default["Mission"]
readmis(x::String) = read_input(x, mis, dmis)
nmisx = readmis("N_missions")
pari = zeros(Int64, iitotal)
parg = zeros(Float64, igtotal)
parm = zeros(Float64, (imtotal, nmisx))
para = zeros(Float64, (iatotal, iptotal, nmisx))
pare = zeros(Float64, (ietotal, iptotal, nmisx))

fuselage = Fuselage()
wing = Wing()
htail = Tail()
vtail = Tail()

# Setup mission variables
ranges = readmis("range")
parm[imRange, :] .= Distance.(ranges)

if !ac_type_fixed
    if parm[imRange, 1] <= 2600 * nmi_to_m
        default = get_default_input_file(0)
        @info """Setting AC Type based on design range: Regional Aircraft"""
    elseif parm[imRange, 1] <= 3115 * nmi_to_m
        default = get_default_input_file(1)
        @info """Setting AC Type based on design range: Narrow Body Aircraft"""
    elseif parm[imRange, 1] <= 8500 * nmi_to_m
        default = get_default_input_file(2)
        @info """Setting AC Type based on design range: Wide Body Aircraft"""
    end    
    ac_type_fixed = true
end

maxpax = readmis("max_payload_in_pax_equivalent") #This represents the maximum aircraft payload in equivalent number of pax
pax = readmis("pax")
exitlimit = readmis("exit_limit")
despax = pax[1] #Design number of passengers

if any(maxpax .< pax)
    error("One or more missions have higher payload than prescribed maximum aircraft payload!:"*
    "\n Max Payload [in pax] = "*string(maxpax)*
    "\n Payloads listed [in pax] = "*string(pax))
end
if any(exitlimit .< pax)
    error("One or more missions have higher passenger counts than the prescribed exit limit!:"*
    "\n Exit limit [in pax] = "*string(exitlimit)*
    "\n Payloads listed [in pax] = "*string(pax))
end

Wpax =  Force(readmis("weight_per_pax"))
parm[imWperpax, :] .= Wpax
parm[imWpay, :] .= pax * Wpax
parg[igWpaymax] = maxpax * Wpax
parg[igfreserve] = readmis("fuel_reserves")
parg[igVne] = Speed(readmis("Vne"))
parg[igNlift] = readmis("Nlift")
fuselage.cabin.design_pax = despax
fuselage.cabin.exit_limit = exitlimit

# Setup option variables
options = read_input("Options", data, default)
doptions = default["Options"]

pari[iiopt] = read_input("optimize", options, doptions)

propsys = read_input("prop_sys_arch", options, doptions)
if lowercase(propsys) == "tf"
    pari[iiengtype] = 1
elseif lowercase(propsys) == "te"
    pari[iiengtype] = 0
else
    error("Propulsion system \"$propsys\" specified. Choose between
    > TF - turbo-fan
    > TE - turbo-electric" )
end

engloc = read_input("engine_location", options, doptions)

if typeof(engloc) == Int
    pari[iiengloc] = engloc
elseif typeof(engloc) <: AbstractString
    engloc = lowercase(engloc)
    if engloc == "wing"
        pari[iiengloc] = 1
    elseif engloc == "fuselage" || engloc == "fuse"
        pari[iiengloc] = 2
    else
        error("Engine location provided is \"$engloc\". Engine position can only be:
        > 1: Engines on \"wing\"
        > 2: Engines on \"fuselage\"")
    end
else
    error("Check engine position input... something isn't right")
end


# Fuel related options
fuel = read_input("Fuel", data, default)
dfuel = default["Fuel"]
readfuel(x::String) = read_input(x, fuel, dfuel)
fueltype = readfuel("fuel_type")
#TODO this needs to be updated once I include Gas.jl into TASOPT
if uppercase(fueltype) == "LH2"
    pari[iifuel] = 40
    
elseif uppercase(fueltype) == "CH4"
    pari[iifuel] = 11
    
elseif uppercase(fueltype) == "JET-A"
    pari[iifuel] = 24

    pare[ieTft, :, :] .= readfuel("fuel_temp") #Temperature of fuel in fuel tank
    pare[ieTfuel, :, :] .= readfuel("fuel_temp") #Initialize fuel temperature as temperature in tank
    parg[igrhofuel] = readfuel("fuel_density")
else
    error("Check fuel type")
end
pari[iifwing]  = readfuel("fuel_in_wing")
pari[iifwcen]  = readfuel("fuel_in_wingcen")
parg[igrWfmax] = readfuel("fuel_usability_factor")
pare[iehvap, :, :] .= readfuel("fuel_enthalpy_vaporization") #Heat of vaporization of the fuel
pare[iehvapcombustor, :, :] .= readfuel("fuel_enthalpy_vaporization") #Heat of vaporization of fuel, if vaporized in combustor

##Takeoff
takeoff = readmis("Takeoff")
dtakeoff = dmis["Takeoff"]
readtakeoff(x) = read_input(x, takeoff, dtakeoff)
parm[imaltTO, :] .= Distance.(readtakeoff("takeoff_alt"))
parg[igmubrake] = readtakeoff("braking_resistance_coeff")
parg[igmuroll]  = readtakeoff("rolling_resistance_coeff")
parg[ighobst]   = Distance(readtakeoff("takeoff_obstacle_height"))
parg[igcdefan]  = readtakeoff("CD_dead_engine")
parg[igCDgear]  = readtakeoff("CD_landing_gear")
parg[igCDspoil] = readtakeoff("CD_spoilers")
parg[iglBFmax]  = Distance(readtakeoff("max_balanced_field_length"))
parg[igNland]   = readtakeoff("Nland")

T0TO = Temp.(readtakeoff("takeoff_T"))
parm[imT0TO, :] .= T0TO 

para[iaclpmax, ipstatic:ipcutback, :] .= readtakeoff("CL_max_perp")
para[iaclpmax, ipclimb1, :] .= readtakeoff("CL_max_perp")
para[iaclpmax, ipdescentn, :] .= readtakeoff("CL_max_perp")

##Climb parameters
climb = readmis("Climb")
dclimb = dmis["Climb"]
parg[iggtocmin] = Angle.(read_input("minimum_top-of-climb_gradient",
                 climb, dclimb))

##Cruise parameters
cruise = readmis("Cruise")
dcruise = dmis["Cruise"]
readcruise(x) = read_input(x, cruise, dcruise)
para[iaalt, ipcruise1, :] .= Distance.(readcruise("cruise_alt"))

para[iaMach, ipclimbn:ipdescent1, :] .= transpose(readcruise("cruise_mach")) #transpose for proper vector broadcasting
para[iaCL, ipclimb1+1:ipdescentn-1, :] .= transpose(readcruise("cruise_CL")) 

##Descent parameters
des = readmis("Descent")
ddes = dmis["Descent"]
readdes(x) = read_input(x, des, ddes)
parm[imgamVDE1, :] .= Angle.(readdes("descent_angle_top-of-descent"))
parm[imgamVDEn, :] .= Angle.(readdes("descent_angle_bottom-of-descent"))

#---------- End Mission vars --------------

# ---------------------------------
# Fuselage
# ---------------------------------
# Setup Fuselage 
fuse = read_input("Fuselage", data, default)
dfuse = default["Fuselage"]

#cabin pressure setting, by explicit pressure value or altitude
#explicit value takes precedence
if  "cabin_pressure" in keys(fuse)
    p_cabin = Pressure.(read_input("cabin_pressure",fuse,dfuse))

else  #if not set explicitly, use altitude (set by default)
    cabinPressureAlt_km = convertDist(parse_unit(read_input("cabin_pressure_altitude",
                                            fuse, dfuse))..., "km")
    _, p_cabin, _, _, _ = atmos(cabinPressureAlt_km)
end
parg[igpcabin] = p_cabin

aero = read_input("Aero", fuse, dfuse)
daero = dfuse["Aero"]
readaero(x) = read_input(x, aero, daero)
    para[iafexcdf, :, :] .= transpose(readaero("excrescence_drag_factor")) #transpose for proper vector broadcasting
    para[iafduo, :, :] .= transpose(readaero("wingroot_fuse_overspeed"))
    para[iafdus, :, :] .= transpose(readaero("wingbreak_fuse_overspeed"))
    para[iafdut, :, :] .= transpose(readaero("wingtip_fuse_overspeed"))

    parg[igCMVf1] = Vol(readaero("fuse_moment_volume_deriv"))
    parg[igCLMf0] = readaero("CL_zero_fuse_moment")
    
    parg[igfBLIf] = readaero("BLI_frac")

weight = read_input("Weights", fuse, dfuse)
dweight = dfuse["Weights"]
readweight(x) = read_input(x, weight, dweight)
    fuselage.weight_frac_frame = readweight("frame")
    fuselage.weight_frac_stringers = readweight("stringer")
    fuselage.weight_frac_skin_addl   = readweight("additional")

    fuselage.fixed.W = Force(readweight("fixed_weight"))

    fuselage.window_W_per_length= readweight("window_per_length")
    fuselage.insulation_W_per_area = readweight("window_insul_per_area")
    fuselage.floor_W_per_area = readweight("floor_weight_per_area")

    fuselage.HPE_sys.W = readweight("HPE_sys_weight_fraction")
    parg[igflgnose] = readweight("LG_nose_weight_fraction")
    parg[igflgmain] = readweight("LG_main_weight_fraction")

    fuselage.APU.W = readweight("APU_weight_fraction")*maxpax*Wpax
    fuselage.seat.W = readweight("seat_weight_fraction")*maxpax*Wpax
    fuselage.added_payload.W = readweight("add_payload_weight_fraction")*maxpax*Wpax

geom = read_input("Geometry", fuse, dfuse)
dgeom = dfuse["Geometry"]
readgeom(x) = read_input(x, geom, dgeom)
    #Boolean to check if cabin length has to be recalculated; if true, this is done 
    #after loading the wing and stabilizer positions
    calculate_cabin = readgeom("calculate_cabin_length") 
    pari[iidoubledeck] = readgeom("double_decker") 

    if pari[iidoubledeck] == 1 #If aircraft is a double decker
        fuselage.cabin.floor_distance = Distance(readgeom("floor_distance")) #read vertical distance between floors
    end

    fuselage.cabin.seat_pitch = Distance(readgeom("seat_pitch"))
    fuselage.cabin.seat_width = Distance(readgeom("seat_width"))
    fuselage.cabin.seat_height = Distance(readgeom("seat_height"))
    fuselage.cabin.aisle_halfwidth = Distance(readgeom("aisle_halfwidth"))
    parg[igrMh] = readgeom("HT_load_fuse_bend_relief")
    parg[igrMv] = readgeom("VT_load_fuse_bend_relief")
    parg[igxlgnose]  = Distance(readgeom("x_nose_landing_gear"))
    parg[igdxlgmain] = Distance(readgeom("x_main_landing_gear_offset"))
    fuselage.APU.r = [Distance(readgeom("x_APU")),0.0,0.0]
    fuselage.HPE_sys.r  = [Distance(readgeom("x_HPE_sys")), 0.0, 0.0]

    fuselage.fixed.r = [Distance(readgeom("x_fixed_weight")),0.0,0.0]

    parg[igxeng] = Distance(readgeom("x_engines"))
    parg[igyeng] = Distance(readgeom("y_critical_engines"))
    
    if readgeom("double_decker")
        fuselage.n_decks =  2
    else
        fuselage.n_decks =  1
    end

    # Number of webs = number of bubbles - 1
    n_bubbles = Int(readgeom("number_of_bubbles"))

    radius = Distance(readgeom("radius"))
    dz = Distance(readgeom("dRadius"))
    dy = Distance(readgeom("y_offset"))

    if n_bubbles > 1 && dy == 0.0
        @warn "Multi-bubble ('$(n_webs+1)') fuselage specified but "*
        "y-offset of bubble set to 0.0. "*
        "Assuming this is a single-bubble design and setting 'number_of_bubbles' = 0"
        n_bubbles = 1
    end

    if n_bubbles == 1
        cross_section = SingleBubble(radius = radius, bubble_lower_downward_shift = dz)
    else
        n_webs = n_bubbles - 1
        cross_section = MultiBubble(radius = radius, bubble_lower_downward_shift = dz,
        bubble_center_y_offset = dy, n_webs = n_webs)
    end

    fuselage.layout.cross_section = cross_section
    fuselage.layout.floor_depth = Distance(readgeom("floor_depth"))
    fuselage.layout.nose_radius = readgeom("a_nose")
    fuselage.layout.tail_radius = readgeom("b_tail")
    fuselage.layout.taper_tailcone = readgeom("tailcone_taper")
    fuse_end = lowercase(readgeom("taper_fuse_to")) 
    if fuse_end == "point"
        fuselage.layout.taper_fuse = 0
    elseif fuse_end == "edge"
        fuselage.layout.taper_fuse = 1
    else
        fuselage.layout.taper_fuse = 0
        @warn "Fuselage can only be closed to a 'point' or an 'edge'"*
                " but '$fuse_end' was provided."*
                " Setting fuselage to end at a point."
    end
    fuselage.layout.x_nose = Distance(readgeom("x_nose_tip")) 
    fuselage.layout.x_pressure_shell_fwd = Distance(readgeom("x_pressure_shell_fwd"))
    fuselage.layout.x_pressure_shell_aft = Distance(readgeom("x_pressure_shell_aft"))
    fuselage.layout.x_start_cylinder = Distance(readgeom("x_start_cylinder"))
    fuselage.layout.x_end_cylinder = Distance(readgeom("x_end_cylinder"))
    fuselage.layout.x_cone_end = Distance(readgeom("x_cone_end"))
    fuselage.layout.x_end = Distance(readgeom("x_end")) 
    fuselage.layout.l_cabin_cylinder = fuselage.layout.x_end_cylinder - fuselage.layout.x_start_cylinder

# ------ End fuse -------


#Fuel storage options
fuse_tank = fuselage_tank() #Initialize struct for fuselage fuel tank params

if pari[iifwing]  == 0 #If fuel is stored in fuselage
    fuel_stor = readfuel("Storage")
    dfuel_stor = dfuel["Storage"]
    readfuel_storage(x::String) = read_input(x, fuel_stor, dfuel_stor)

    fuse_tank.placement = readfuel_storage("tank_placement")
    fuse_tank.fueltype = fueltype
    fuse_tank.clearance_fuse = Distance(readfuel_storage("fuselage_clearance"))

    fuse_tank.size_insulation = readfuel_storage("size_insulation")
    fuse_tank.t_insul = readfuel_storage("insulation_segment_base_thickness")
    insul_mats_names = readfuel_storage("insulation_material")
    insul_mats = []
    for insul_mat_name in insul_mats_names
        push!(insul_mats, ThermalInsulator(insul_mat_name))
    end
    fuse_tank.material_insul = insul_mats
    if fuse_tank.size_insulation
        fuse_tank.boiloff_rate = readfuel_storage("cruise_boiloff_rate")
        fuse_tank.iinsuldes = readfuel_storage("insulation_thicknesses_design_indices")
    end
    
    inner_mat_name = readfuel_storage("inner_vessel_material")
    fuse_tank.inner_material = StructuralAlloy(inner_mat_name)
    
    fuse_tank.ARtank = readfuel_storage("tank_aspect_ratio")
    fuse_tank.theta_inner = Angle(readfuel_storage("inner_vessel_support_angle"))

    fuse_tank.pvent = Pressure(readfuel_storage("pressure_venting"))
    fuse_tank.pinitial = Pressure(readfuel_storage("pressure_initial"))
    fuse_tank.t_hold_orig = Time(readfuel_storage("hold_departure"))
    fuse_tank.t_hold_dest = Time(readfuel_storage("hold_arrival"))
    
    fuse_tank.ftankadd = readfuel_storage("additional_mass_fraction")
    fuse_tank.ew = readfuel_storage("weld_efficiency")
    fuse_tank.ullage_frac = readfuel_storage("ullage_fraction")
    fuse_tank.qfac = readfuel_storage("heat_leak_factor")
    fuse_tank.TSLtank = Temp(readfuel_storage("SL_temperature_for_tank"))
    fuse_tank.pfac = readfuel_storage("pressure_rise_factor")

    flag_vacuum = TASOPT.CryoTank.check_vacuum(fuse_tank.material_insul) #flag to check if an outer vessel is needed

    if flag_vacuum #If tank is double-walled
        outer_mat_name = readfuel_storage("outer_vessel_material")
        fuse_tank.outer_material = StructuralAlloy(outer_mat_name)

        theta_outer_str = readfuel_storage("outer_vessel_support_angles")
        theta_outer = []
        for θstr in theta_outer_str
            push!(theta_outer,  Angle(θstr))
        end
        fuse_tank.theta_outer = theta_outer
        fuse_tank.Ninterm = 1.0 #Initial guess for first iteration
    end

    #Find number of tanks from placement
    if (fuse_tank.placement == "front") || (fuse_tank.placement == "rear")
        pari[iinftanks] = 1
    elseif (fuse_tank.placement == "both") 
        pari[iinftanks] = 2
    end
end
# ---------------------------------
# Wing
# ---------------------------------
# Setup wing
wing_i = read_input("Wing", data, default)
dwing = default["Wing"]
readwing(x) = read_input(x, wing_i, dwing)
    wing.planform = readwing("wing_planform")
    wing.has_strut = readwing("strut_braced_wing")

    wing.layout.sweep = readwing("sweep")
    wing.layout.AR = readwing("AR")
    wing.layout.max_span = Distance(readwing("maxSpan"))

    wing.inboard.λ = readwing("inner_panel_taper_ratio")
    wing.outboard.λ = readwing("outer_panel_taper_ratio")
    wing.layout.ηs    = readwing("panel_break_location")
    if !(0 ≤ wing.layout.ηs ≤ 1.0)
        @warn "Wing span break location input was $(wing.layout.ηs); ηs must be 0 ≤ ηs ≤ 1.0"
        if wing.layout.ηs > 1.0
            wing.layout.ηs = 1.0
        else 
            wing.layout.ηs = 0.0
        end
        @warn "ηs set to $(wing.layout.ηs)"
    end

    wing.layout.root_span = 2*Distance(readwing("center_box_halfspan"))
    wing.inboard.cross_section.width_to_chord  = readwing("box_width_to_chord")
    wing.outboard.cross_section.width_to_chord  = readwing("box_width_to_chord")
   
    wing.inboard.cross_section.thickness_to_chord = readwing("root_thickness_to_chord")
    wing.outboard.cross_section.thickness_to_chord = readwing("spanbreak_thickness_to_chord")
    
    wing.inboard.cross_section.web_to_box_height  = readwing("hweb_to_hbox")
    wing.outboard.cross_section.web_to_box_height = readwing("hweb_to_hbox")
    wing.layout.spar_box_x_c = readwing("spar_box_x_c")

    wing.layout.box_x = Distance(readwing("x_wing_box"))
    wing.layout.z = Distance(readwing("z_wing"))

    parg[igdxeng2wbox] = wing.layout.box_x - parg[igxeng] #TODO add this as a function of wing

    ## Strut details only used if strut_braced_wing is true
    if wing.has_strut
        wing.strut.z  = Distance(readwing("z_strut"))
        wing.strut.thickness_to_chord  = readwing("strut_toc")
        wing.strut.local_velocity_ratio = readwing("strut_local_velocity_ratio")
    end

    airfoil_data = joinpath(__TASOPTroot__,"airfoil_data/", readwing("airfoil"))
    wing.airsection = TASOPT.aerodynamics.airtable(airfoil_data);

# ----------------------------------
# ------- Wing Aerodynamics --------
# ----------------------------------
aero = readwing("Aero")
daero = dwing["Aero"]
    wing.fuse_lift_carryover = readaero("fuselage_lift_carryover_loss_factor")
    wing.tip_lift_loss = readaero("wing_tip_lift_rolloff_factor")

    para[iacdfw, 1:iptotal, :]   .= readaero("lowspeed_cdf")  #  cdfw    wing profile cd for low speed (takeoff, initial climb)
    para[iacdpw, 1:iptotal, :]   .= readaero("lowspeed_cdp")  #  cdpw    
    para[iaRerefw, 1:iptotal, :] .= readaero("Re_ref")  #  Rerefw

    para[iacdfs, 1:iptotal, :]   .= readaero("strut_lowspeed_cdf")  #  cdfs    strut profile cd (not used if there's no strut)
    para[iacdps, 1:iptotal, :]   .= readaero("strut_lowspeed_cdp")  #  cdps    
    para[iaRerefs, 1:iptotal, :] .= readaero("strut_Re_ref")        #  Rerefs  

    para[iaaRexp, 1:iptotal, :] .= readaero("Reynolds_scaling")
    para[iafexcdw, 1:iptotal, :] .= readaero("excrescence_drag_factor")
    parg[igfBLIw] = readaero("BLI_frac")

#- wing spanwise cl and cm distributions over mission
#- ( rclo = clo/clo = 1.0  by definition, so it's not specified )
#- takeoff, initial climb
takeoff = readaero("Takeoff")
dtakeoff = daero["Takeoff"]
    #transpose for proper vector broadcasting
    para[iarcls, 1:ipclimb1, :] .= transpose(readtakeoff("cls_clo"))    #  rcls    break/root cl ratio = cls/clo
    para[iarclt, 1:ipclimb1, :] .= transpose(readtakeoff("clt_clo"))    #  rclt    tip  /root cl ratio = clt/clo
    para[iacmpo, 1:ipclimb1, :] .= transpose(readtakeoff("cm_o"))      #  cmpo    root  cm
    para[iacmps, 1:ipclimb1, :] .= transpose(readtakeoff("cm_s"))      #  cmps    break cm
    para[iacmpt, 1:ipclimb1, :] .= transpose(readtakeoff("cm_t"))      #  cmpt    tip   cm

# Clean climb cruise descent and for wing structure sizing
climb = readaero("Climb")
dclimb = daero["Climb"]
readclimb(x) = read_input(x, climb, dclimb)
    #transpose for proper vector broadcasting
    para[iarcls, ipclimb1+1:ipdescentn-1, :] .= transpose(readclimb("cls_clo"))   #  rcls    break/root cl ratio = cls/clo
    para[iarclt, ipclimb1+1:ipdescentn-1, :] .= transpose(readclimb("clt_clo"))   #  rclt    tip  /root cl ratio = clt/clo
    para[iacmpo, ipclimb1+1:ipdescentn-1, :] .= transpose(readclimb("cm_o"))      #  cmpo    root  cm
    para[iacmps, ipclimb1+1:ipdescentn-1, :] .= transpose(readclimb("cm_s"))      #  cmps    break cm
    para[iacmpt, ipclimb1+1:ipdescentn-1, :] .= transpose(readclimb("cm_t"))      #  cmpt    tip   cm

# Landing, forward CG tail sizing case
land = readaero("Landing")
dland = daero["Landing"]
readland(x) = read_input(x, land, dland)
    #transpose for proper vector broadcasting
    para[iarcls, ipdescentn, :] .= transpose(readland("cls_clo"))   #  rcls    break/root cl ratio = cls/clo
    para[iarclt, ipdescentn, :] .= transpose(readland("clt_clo"))   #  rclt    tip  /root cl ratio = clt/clo
    para[iacmpo, ipdescentn, :] .= transpose(readland("cm_o"))      #  cmpo    root  cm
    para[iacmps, ipdescentn, :] .= transpose(readland("cm_s"))      #  cmps    break cm
    para[iacmpt, ipdescentn, :] .= transpose(readland("cm_t"))      #  cmpt    tip   cm
    
# ----------------------------------
# ---------- Wing Weight -----------
# ----------------------------------
# Wing weight fractions of flight surfaces
# and secondary wing components, 
weight = readwing("Weightfracs")
dweight = dwing["Weightfracs"]
    wing.weight_frac_flap  = readweight("flap")
    wing.weight_frac_slat = readweight("slat")
    wing.weight_frac_ailerons = readweight("aileron")
    wing.weight_frac_leading_trailing_edge = readweight("leading_trailing_edge")
    wing.weight_frac_ribs = readweight("ribs")
    wing.weight_frac_spoilers = readweight("spoilers")
    wing.weight_frac_attachments = readweight("attachments")

# ---- End Wing -----

# ---------------------------------
# Stabilizers
# ---------------------------------
tails = read_input("Stabilizers", data, default)
dtails = default["Stabilizers"]
readtails(x) = read_input(x, tails, dtails)
    #transpose for proper broadcasting
    para[iacdft, 1:iptotal, :]   .= transpose(readtails("lowspeed_cdf"))  #  cdft    tail profile cd
    para[iacdpt, 1:iptotal, :]   .= transpose(readtails("lowspeed_cdp"))  #  cdpt    
    para[iaRereft, 1:iptotal, :] .= transpose(readtails("Re_ref"))  #  Rereft  

    para[iafexcdt, 1:iptotal, :] .= transpose(readtails("excrescence_drag_factor"))

    htail_input = readtails("Htail")
    dhtail = dtails["Htail"]

readhtail(x) = read_input(x, htail_input, dhtail)
    htail.layout.AR = readhtail("AR_Htail")
    htail.outboard.λ = readhtail("taper")
    htail.layout.sweep = readhtail("sweep")
    htail.layout.root_span = 2*Distance(readhtail("center_box_halfspan"))

    htail.layout.box_x  = Distance(readhtail("x_Htail"))
    htail.layout.z = Distance(readhtail("z_Htail"))
    htail.ntails  = readhtail("number_Htails")

    htail.CL_CLmax = readhtail("max_tail_download")

    htail_size = lowercase(readhtail("HTsize"))
    if htail_size == "vh"
        htail.size = 1
        htail.volume = readhtail("Vh")
    elseif htail_size == "maxforwardcg"
        htail.size = 2
        htail.CL_max_fwd_CG = readhtail("CLh_at_max_forward_CG")
        htail.volume = 1.0
    else
        error("Horizontal tail can only be sized via:
            1: specified tail volume coeff \"Vh\";
            2: specified CLh at max-forward CG case during landing (\"maxforwardCG\")")
    end


    movewing = readhtail("move_wingbox")
    if typeof(movewing) == Int
        htail.move_wingbox = movewing
    elseif typeof(movewing) <: AbstractString
        movewing = lowercase(movewing)
        if movewing =="fix"
            htail.move_wingbox = 0
        elseif movewing == "clhspec"
            htail.move_wingbox = 1
        elseif movewing == "smmin"
            htail.move_wingbox = 2
        else
            error("Wing position during horizontal tail sizing can only be sized via:
            0: \"fix\" wing position;
            1: move wing to get CLh=\"CLhspec\" in cruise 
            2: move wing to get min static margin = \"SMmin\"")
        end
    else
        error("Check wing position input during htail sizing... something isn't right")
    end
    htail.SM_min = readhtail("SM_min")

    parg[igCLhspec] = readhtail("CLh_spec")

    htail.downwash_factor = readhtail("downwash_factor")
    parg[igdCLnda] = readhtail("nacelle_lift_curve_slope")

    parg[igfCDhcen] = readhtail("CD_Htail_from_center")
    htail.CL_max  = readhtail("CLh_max")

    htail.weight_fraction_added = readhtail("added_weight_fraction")

    htail.outboard.cross_section.width_to_chord = readhtail("box_width_to_chord")
    htail.outboard.cross_section.thickness_to_chord = readhtail("box_height_chord")
    htail.outboard.cross_section.web_to_box_height  = readhtail("web_height_hbox")


vtail_input = readtails("Vtail")
dvtail = dtails["Vtail"]
readvtail(x) = read_input(x, vtail_input, dvtail)
    vtail.layout.AR = readvtail("AR_Vtail")
    vtail.outboard.λ = readvtail("taper")
    vtail.layout.sweep  = readvtail("sweep")
    vtail.layout.root_span = Distance(readvtail("center_box_halfspan"))
    vtail.layout.box_x  = Distance(readvtail("x_Vtail"))
    vtail.ntails  = readvtail("number_Vtails")

    vtail_size = lowercase(readvtail("VTsize"))
    if vtail_size == "vv"
        vtail.size = 1
        vtail.volume = readvtail("Vv")
    elseif vtail_size == "oei"
        vtail.size = 2
        parg[igCLveout] = readvtail("CLv_at_engine_out")
    else
        error("Vertical tail can only be sized via:
            1: specified tail volume coeff \"Vv\";
            2: specified CL at one engine out trim (\"OEI\")")
    end

    vtail.CL_max = readvtail("CLv_max")

    vtail.weight_fraction_added = readvtail("added_weight_fraction")
    vtail.outboard.cross_section.width_to_chord = readvtail("box_width_to_chord")
    vtail.outboard.cross_section.thickness_to_chord = readvtail("box_height_chord")
    vtail.outboard.cross_section.web_to_box_height  = readvtail("web_height_hbox")

# ----- End Stabilizers -----

# ---------------------------------
# Recalculate cabin length
if calculate_cabin #Resize the cabin if desired, keeping deltas
    @info "Fuselage and stabilizer layouts have been overwritten; deltas will be maintained."
    update_fuse_for_pax!(pari, parg, fuselage, fuse_tank, wing, htail, vtail) #update fuselage dimensions
end
# ---------------------------------

# ---------------------------------
# Structures
# ---------------------------------

structures = read_input("Structures", data, default)
dstructures = default["Structures"]
readstruct(x) = read_input(x, structures, dstructures)
    parg[igsigfac] = readstruct("stress_factor")

    #- fuselage shell modulus ratio, for bending material sizing
    fuselage.ratio_young_mod_fuse_bending = readstruct("fuse_shell_modulus_ratio")

    caps_max_avg = readstruct("caps_max_avg_stress")
    caps_safety_fac = readstruct("caps_safety_factor")
    wing.inboard.caps.material = StructuralAlloy(readstruct("caps_material"),
                                max_avg_stress = caps_max_avg,
                                safety_factor = caps_safety_fac)
    wing.outboard.caps.material = StructuralAlloy(readstruct("caps_material"),
                                max_avg_stress = caps_max_avg,
                                safety_factor = caps_safety_fac)

    webs_max_avg = readstruct("webs_max_avg_stress")
    webs_safety_fac = readstruct("webs_safety_factor")
    wing.inboard.webs.material = StructuralAlloy(readstruct("webs_material"),
                                max_avg_stress = webs_max_avg,
                                safety_factor = webs_safety_fac)
    wing.outboard.webs.material = StructuralAlloy(readstruct("webs_material"),
                                max_avg_stress = webs_max_avg,
                                safety_factor = webs_safety_fac)

    skin_max_avg = readstruct("skin_max_avg_stress")
    skin_safety_fac = readstruct("skin_safety_factor")
    fuselage.skin.material =  StructuralAlloy(readstruct("skin_material"), 
                                    max_avg_stress = skin_max_avg,
                                    safety_factor = skin_safety_fac)

    cone_max_avg = readstruct("cone_max_avg_stress")
    cone_safety_fac = readstruct("cone_safety_factor")
    fuselage.cone.material =  StructuralAlloy(readstruct("cone_material"), 
                                    max_avg_stress = cone_max_avg,
                                    safety_factor = cone_safety_fac)

    bend_max_avg = readstruct("bending_max_avg_stress")
    bend_safety_fac = readstruct("bending_safety_factor")
    fuselage.bendingmaterial_h.material = StructuralAlloy(readstruct("bending_material"),
                                max_avg_stress = bend_max_avg,
                                safety_factor = bend_safety_fac)

    fuselage.bendingmaterial_v.material = fuselage.bendingmaterial_h.material                            

    floor_max_avg = readstruct("floor_max_avg_stress")
    floor_safety_fac = readstruct("floor_safety_factor")
    fuselage.floor.material = StructuralAlloy(readstruct("floor_material"),
        max_avg_stress=floor_max_avg,
        safety_factor=floor_safety_fac)

# ---------------------------------
# Propulsion systems
# ---------------------------------

prop = read_input("Propulsion", data, default)
dprop = default["Propulsion"]
readprop(x) = read_input(x, prop, dprop)
    parg[igneng] = readprop("number_of_engines")
    parg[igTmetal] = Temp.(readprop("T_max_metal"))
    parg[igfTt4CL1] = readprop("Tt4_frac_bottom_of_climb")
    parg[igfTt4CLn] = readprop("Tt4_frac_top_of_climb")

    pare[ieTt4, :, :] .= transpose(Temp.(readprop("Tt4_cruise"))) #transpose for proper broadcasting

    Tt4TO = transpose(Temp.(readprop("Tt4_takeoff")))
    pare[ieTt4, ipstatic, :] .= Tt4TO
    pare[ieTt4, iprotate, :] .= Tt4TO
    pare[ieTt4, iptakeoff, :] .= Tt4TO

    pare[ieT0, ipstatic, :] .= T0TO
    pare[ieT0, iprotate, :] .= T0TO
    pare[ieT0, iptakeoff, :] .= T0TO

    # Core in clean-flow -> 0; Core ingests KE defect -> 1
    pari[iiBLIc] = !readprop("core_in_clean_flow")

    #Turbomachinery
    turb = readprop("Turbomachinery")
    dturb = dprop["Turbomachinery"]

readturb(x) = read_input(x, turb, dturb)
    Gearf = readturb("gear_ratio")
    BPR = readturb("BPR")
    OPR = readturb("OPR")
    pif = readturb("Fan_PR")
    pilc = readturb("LPC_PR")
    pihc = OPR./pilc

    pid = readturb("diffuser_PR")
    pib = readturb("burner_PR")
    pifn = readturb("fan_nozzle_PR")
    pitn = readturb("core_nozzle_PR")

    epolf  = readturb("fan_eta_poly") 
    epollc = readturb("LPC_eta_poly") 
    epolhc = readturb("HPC_eta_poly") 
    epolht = readturb("HPT_eta_poly") 
    epollt = readturb("LPT_eta_poly") 

    pifK = readturb("FPR0")
    epfK = readturb("Kf_polyeff")
    HTRf  = readturb("HTR_fan")
    HTRlc = readturb("HTR_LPC")
    HTRhc = readturb("HTR_HPC")

    M2  = readturb("M2")
    M25 = readturb("M25")

    epsl = readturb("low_spool_loss")
    epsh = readturb("high_spool_loss")

comb = read_input("Combustor", prop, dprop)
dcomb = dprop["Combustor"]
    etab = read_input("combustion_efficiency", comb, dcomb)

pare[iepid, :, :] .= pid
pare[iepib, :, :] .= pib
pare[iepifn, :, :] .= pifn
pare[iepitn, :, :] .= pitn
pare[iepif, :, :] .= pif
pare[iepilc, :, :] .= pilc
pare[iepihc, :, :] .= pihc
pare[ieepolf, :, :] .= epolf
pare[ieepollc, :, :] .= epollc
pare[ieepolhc, :, :] .= epolhc
pare[ieepolht, :, :] .= epolht
pare[ieepollt, :, :] .= epollt
pare[ieetab, :, :] .= etab
pare[iepifK, :, :] .= pifK
pare[ieepfK, :, :] .= epfK
pare[ieBPR, :, :] .= BPR
pare[ieM2, :, :] .= M2
pare[ieM25, :, :] .= M25
pare[ieepsl, :, :] .= epsl
pare[ieepsh, :, :] .= epsh

parg[igGearf] = Gearf
parg[igHTRf] = HTRf
parg[igHTRlc] = HTRlc
parg[igHTRhc] = HTRhc

# Cooling
cool = readprop("Cooling")
dcool = dprop["Cooling"]
readcool(x) = read_input(x, cool, dcool)
    dTstrk = Temp(readcool("hot_streak_T_allowance"))
    Mtexit = readcool("M_turbine_blade_exit")
    StA = readcool("St")

    efilm = readcool("e_film_cooling")
    tfilm = readcool("t_film_cooling")

    M41 = readcool("M41")
    ruc = readcool("cooling_air_V_ratio")

    pare[ieM4a, :, :] .= M41
    pare[ieruc, :, :] .= ruc
    pare[iedTstrk, :, :] .= dTstrk
    pare[ieMtexit, :, :] .= Mtexit
    pare[ieStA, :, :] .= StA
    pare[ieefilm, :, :] .= efilm
    pare[ietfilm, :, :] .= tfilm

# Offtakes
off = readprop("Offtakes")
doff = dprop["Offtakes"]
readoff(x) = read_input(x, off, doff)
    mofftpax  = readoff("LPC_mass_offtake_per_pax")
    mofftmMTO = readoff("LPC_mass_offtake_per_max_mass")

    Pofftpax  = readoff("Low_spool_power_offtake_per_pax")
    PofftmMTO = readoff("Low_spool_power_offtake_per_max_mass")

    Ttdischarge = readoff("Tt_offtake_air")
    Ptdischarge = readoff("Pt_offtake_air")

    # TODO Tt9 is really a terrible numbering convention for the discharge temp 

    pare[ieTt9, :, :] .= Ttdischarge
    pare[iept9, :, :] .= Ptdischarge

    parg[igmofWpay] = mofftpax ./ parm[imWperpax, 1]
    parg[igmofWMTO] = mofftmMTO / gee
    parg[igPofWpay] = Pofftpax ./ parm[imWperpax, 1]
    parg[igPofWMTO] = PofftmMTO / gee

## Nozzle areas
noz = readprop("Nozzles")
dnoz = dprop["Nozzles"]
corenoz = read_input("core_nozzle_area", noz, dnoz)
dcorenoz = dnoz["core_nozzle_area"]
readcnoz(x) = read_input(x, corenoz, dcorenoz)
    A5static   = readcnoz("static")
    A5takeoff  = readcnoz("rotation")
    A5cutback  = readcnoz("cutback")
    A5climb1   = readcnoz("climbstart")
    A5climbn   = readcnoz("climbend")
    A5descent1 = readcnoz("descentstart")
    A5descentn = readcnoz("descentend")

fannoz = read_input("fan_nozzle_area", noz, dnoz)
dfannoz = dnoz["fan_nozzle_area"]
readfnoz(x) = read_input(x, fannoz, dfannoz)
    A7static   = readfnoz("static")
    A7takeoff  = readfnoz("rotation")
    A7cutback  = readfnoz("cutback")
    A7climb1   = readfnoz("climbstart")
    A7climbn   = readfnoz("climbend")
    A7descent1 = readfnoz("descentstart")
    A7descentn = readfnoz("descentend")

    pare[ieA7fac, ipstatic, :] .= A7static
    pare[ieA7fac, iprotate, :] .= A7takeoff
    pare[ieA7fac, iptakeoff, :] .= A7takeoff
    pare[ieA7fac, ipcutback, :] .= A7cutback

    pare[ieA5fac, ipstatic, :] .= A5static
    pare[ieA5fac, iprotate, :] .= A5takeoff
    pare[ieA5fac, iptakeoff, :] .= A5takeoff
    pare[ieA5fac, ipcutback, :] .= A5cutback

    for ip = ipclimb1:ipclimbn

        frac = (ip - ipclimb1) /  (ipclimbn - ipclimb1)

        pare[ieA7fac, ip, :] .= A7climb1 * (1.0 - frac) + A7climbn * frac
        pare[ieA5fac, ip, :] .= A5climb1 * (1.0 - frac) + A5climbn * frac

    end

    pare[ieA7fac, ipcruise1:ipcruisen, :] .= 1.0
    pare[ieA5fac, ipcruise1:ipcruisen, :] .= 1.0

    for ip = ipdescent1:ipdescentn

        frac = (ip - ipdescent1) / (ipdescentn - ipdescent1)

        pare[ieA7fac, ip, :] .= A7descent1 * (1.0 - frac) + A7descentn * frac
        pare[ieA5fac, ip, :] .= A5descent1 * (1.0 - frac) + A5descentn * frac

    end

    pare[ieA7fac, iptest, :] .= A7static
    pare[ieA5fac, iptest, :] .= A5static

nac = readprop("Nacelles")
dnac = dprop["Nacelles"]
    #- nacelle drag stuff
    parg[igrSnace] = read_input("nacelle_pylon_wetted_area_ratio", nac, dnac)
    parg[igrVnace] = read_input("nacelle_local_velocity_ratio", nac, dnac)

weight = readprop("Weight")
dweight = dprop["Weight"]
    parg[igfeadd] = read_input("engine_access_weight_fraction", weight, dweight)
    parg[igfpylon] = read_input("pylon_weight_fraction", weight, dweight)
    TF_wmodel = read_input("weight_model", weight, dweight)
    if lowercase(TF_wmodel) == "md"
        pari[iiengwgt] = 0
    elseif lowercase(TF_wmodel) == "basic"
        pari[iiengwgt] = 1
    elseif lowercase(TF_wmodel) == "advanced"
        pari[iiengwgt] = 2
    else
        error("\"$TF_wmodel\" engine weight model was specifed. 
        Engine weight can only be \"MD\", \"basic\" or \"advanced\".")
    end

# Heat exchangers

HEx = readprop("HeatExchangers")
dHEx = dprop["HeatExchangers"]
    parg[igHXaddmassfrac] = read_input("added_mass_frac", HEx, dHEx)

    pare[iefrecirc, :, :] .= read_input("recirculation_flag", HEx, dHEx)
    pare[ierecircT, :, :] .= read_input("recirculation_temperature", HEx, dHEx)
    pare[ieDi, :, :] .= read_input("core_inner_diameter", HEx, dHEx)
    
    pare[iePreCorder, :, :] .= read_input("precooler_order", HEx, dHEx)
    pare[iePreCepsilon, :, :] .= read_input("precooler_effectiveness", HEx, dHEx)
    pare[iePreCMp, :, :] .= read_input("precooler_inlet_mach", HEx, dHEx)

    pare[ieInterCorder, :, :] .= read_input("intercooler_order", HEx, dHEx)
    pare[ieInterCepsilon, :, :] .= read_input("intercooler_effectiveness", HEx, dHEx)
    pare[ieInterCMp, :, :] .= read_input("intercooler_inlet_mach", HEx, dHEx)

    pare[ieRegenorder, :, :] .= read_input("regenerative_order", HEx, dHEx)
    pare[ieRegenepsilon, :, :] .= read_input("regenerative_effectiveness", HEx, dHEx)
    pare[ieRegenMp, :, :] .= read_input("regenerative_inlet_mach", HEx, dHEx)

    pare[ieTurbCorder, :, :] .= read_input("turbine_cooler_order", HEx, dHEx)
    pare[ieTurbCepsilon, :, :] .= read_input("turbine_cooler_effectiveness", HEx, dHEx)
    pare[ieTurbCMp, :, :] .= read_input("turbine_cooler_inlet_mach", HEx, dHEx)


return TASOPT.aircraft(name, description,
    pari, parg, parm, para, pare, [false], fuse_tank, fuselage, wing, htail, vtail)

end

function load_default_model()
    println("Loading default aircraft model")
    read_aircraft_model()
end