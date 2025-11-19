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

function get_template_input_file(designrange)
    if designrange <= 2600 * nmi_to_m
        templatefile = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml")
    elseif designrange <= 3115 * nmi_to_m
        templatefile = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml")
    elseif designrange <= 8500 * nmi_to_m
        templatefile = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_wide.toml")
    else
        println("\n")
        @warn """Requested aircraft design range exceeds expected capability. Selecting Wide Body Aircraft Template, but be warned. """
        templatefile = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_wide.toml")
    end
    return templatefile
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
Times(x)     = convertTime(parse_unit(x)...)
Temp(x)     = convertTemp(parse_unit(x)...)


"""
    read_aircraft_model(datafile; 
    defaultfile = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml"))

Reads a specified TOML file that describes a TASOPT `aircraft` model 
with a fall back to the default `aircraft` definition 
provided in \"/example/defaults/default_input.toml\""

!!! note "Deviating from default"
    Extending `read_input.jl` and `save_model.jl` is recommended for models deviating appreciably 
    from the default functionality. Thorough knowledge of the model is required.

# Examples
```julia-repl
julia> read_aircraft_model("examples/defaults/default_input.toml")


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
    datafile=joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml"); 
    templatefile = "")

data = TOML.parsefile(datafile)

# Get template input file, with appropriate user notices when needed
# handle default templatefile value
if templatefile == ""
    @info "No template input file provided. Proceeding with template file as determined by design mission range."
    templatefile = nothing
# check if provided template input file is extant
elseif !isfile(templatefile)
    #if not, warn that we're ignoring it
    @warn "Template input file provided does not exist: $templatefile \n Proceeding with template file as determined by design mission range."
    templatefile = nothing
end
#if no valid templatefile provided
if isnothing(templatefile)
    #determine appropriate template input file based on mission range
    designrange = Distance.(data["Mission"]["range"])[1] #in meters
    templatefile = get_template_input_file(designrange)
    @info "Template input file selected: $templatefile"
end

default = TOML.parsefile(templatefile)
ac_descrip = get(data, "Aircraft Description", Dict{})
name = get(ac_descrip, "name", "Untitled Model")
description = get(ac_descrip, "description", "---")
is_sized = get(ac_descrip, "is_sized",[false])


#Get number of missions to create data arrays
mis = read_input("Mission", data, default)
dmis = default["Mission"]
readmis(x::String) = read_input(x, mis, dmis)
nmisx = readmis("N_missions")
parg = zeros(Float64, igtotal)
parm = zeros(Float64, (imtotal, nmisx))
para = zeros(Float64, (iatotal, iptotal, nmisx))
pare = zeros(Float64, (ietotal, iptotal, nmisx))

wing = Wing()
htail = Tail()
vtail = Tail()
landing_gear = LandingGear()

# Setup mission variables
ranges = readmis("range")
parm[imRange, :] .= Distance.(ranges)
Wpax =  Force(readmis("weight_per_pax"))

#Weight() can take in "pax" as a unit; the weight of a passenger is user-defined
function Weight(x)
    if x isa AbstractVector
        return [Weight(w) for w in x]  # Recursively call Weight on each element
    elseif x isa AbstractString
        value, unit = parse_unit(x)
        if unit == "pax" 
            return value * Wpax
        else 
            return convertForce(value, unit)
        end
    elseif x isa Float64
        return x
    else
        throw(ArgumentError("Unsupported input type: $(typeof(x))"))
    end
end

maxpay = Weight(readmis("max_payload")) #This represents the maximum aircraft payload
                                #This may exceed the seatable capacity to account for belly cargo

payload = Weight(readmis("payload"))
exitlimit = readmis("exit_limit") #Maximum number of pax that could fit in cabin in an all-economy layout

parm[imWperpax, :] .= Wpax
parm[imWpay, :] .= payload
parg[igWpaymax] = maxpay 
parg[igfreserve] = readmis("fuel_reserves")
parg[igVne] = Speed(readmis("Vne"))
parg[igNlift] = readmis("Nlift")

# Setup option variables
options = read_input("Options", data, default)
doptions = default["Options"]


# -----------------------------
# Engine model setup
# ------------------------------
engloc = read_input("engine_location", options, doptions)
calculate_takeoff = true #true by default

#throw error if engloc isn't a string indicating a supported location
if !(typeof(engloc) <: AbstractString && engloc in ["wing", "fuselage"])
   error("Engine location provided is \"$engloc\". Engine position can only be:
        > \"wing\" - engines under wing
        > \"fuselage\" - engines on aft fuselage")
end

trefftz_resolution_str = read_input("trefftz_resolution", options, doptions)
trefftz_k_tip = read_input("trefftz_k_tip", options, doptions)
trefftz_bunch = read_input("trefftz_bunch", options, doptions)
wing_root_contraction = read_input("wing_root_contraction", options, doptions)
tail_root_contraction = read_input("tail_root_contraction", options, doptions)
TREFFTZ_CONFIG = aerodynamics.get_trefftz_config(trefftz_resolution_str;
                                                  k_tip=trefftz_k_tip,
                                                  bunch=trefftz_bunch,
                                                  wing_root_contraction=wing_root_contraction,
                                                  tail_root_contraction=tail_root_contraction)


# Fuel related options
fuel = read_input("Fuel", data, default)
dfuel = default["Fuel"]
readfuel(x::String) = read_input(x, fuel, dfuel)
fueltype = readfuel("fuel_type")
#TODO this needs to be updated once Prash includes Gas.jl into TASOPT

#check input, perform actions based on fuel type
if compare_strings(fueltype, "JET-A")
    pare[ieTft, :, :] .= readfuel("fuel_temp") #Temperature of fuel in fuel tank
    pare[ieTfuel, :, :] .= readfuel("fuel_temp") #Initialize fuel temperature as temperature in tank
    parg[igrhofuel] = readfuel("fuel_density")
    ifuel = 24
elseif compare_strings(fueltype, "LH2") 
    ifuel = 40
elseif compare_strings(fueltype, "CH4")
    ifuel = 11
#throw error if fueltype isn't a supported value
else 
    error("'$fueltype' is not a supported fuel type (e.g., \"JET-A\", \"LH2\", \"CH4\")")
end

has_centerbox_fuel  = readfuel("fuel_in_wingcen")
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

geom = read_input("Geometry", fuse, dfuse)
dgeom = dfuse["Geometry"]
readgeom(x) = read_input(x, geom, dgeom)
    
    # Before doing anything else, first learn the shape of the fuselage 
    # cross-section so you can create the appropriate typed Fuselage.
    
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
    fuselage = Fuselage{typeof(cross_section)}() #Create the right type of fuselage
    
    fuselage.layout.cross_section = cross_section
    fuselage.cabin.exit_limit = exitlimit
    #Boolean to check if cabin length has to be recalculated; if true, this is done 
    #after loading the wing and stabilizer positions
    calculate_cabin = readgeom("calculate_cabin_length") 
    is_doubledecker = Bool(readgeom("double_decker"))

    if is_doubledecker #If aircraft is a double decker
        fuselage.n_decks =  2
        fuselage.cabin.floor_distance = Distance(readgeom("floor_distance")) #read vertical distance between floors
        fuselage.cabin.unit_load_device = readgeom("unit_load_device")
        fuselage.cabin.min_top_cabin_height = Distance(readgeom("min_top_cabin_height"))
    else
        fuselage.n_decks =  1
    end
    if calculate_cabin
        fuselage.cabin.front_seat_offset = Distance(readgeom("front_seat_offset"))
        fuselage.cabin.rear_seat_offset = Distance(readgeom("rear_seat_offset"))
    end

    fuselage.cabin.seat_pitch = Distance(readgeom("seat_pitch"))
    fuselage.cabin.seat_width = Distance(readgeom("seat_width"))
    fuselage.cabin.seat_height = Distance(readgeom("seat_height"))
    fuselage.cabin.aisle_halfwidth = Distance(readgeom("aisle_halfwidth"))
    parg[igrMh] = readgeom("HT_load_fuse_bend_relief")
    parg[igrMv] = readgeom("VT_load_fuse_bend_relief")
    fuselage.APU.r = [Distance(readgeom("x_APU")),0.0,0.0]
    fuselage.HPE_sys.r  = [Distance(readgeom("x_HPE_sys")), 0.0, 0.0]

    fuselage.fixed.r = [Distance(readgeom("x_fixed_weight")),0.0,0.0]

    parg[igxeng] = Distance(readgeom("x_engines"))
    parg[igyeng] = Distance(readgeom("y_critical_engines"))

    
    fuselage.layout.floor_depth = Distance(readgeom("floor_depth"))
    fuselage.layout.nose_radius = readgeom("a_nose")
    fuselage.layout.tail_radius = readgeom("b_tail")
    fuselage.layout.taper_tailcone = readgeom("tailcone_taper")

    fuse_end = readgeom("tapers_to")
    #throw error if fuse_end isn't a supported fuse taper
    if !(fuse_end in ["point", "edge"])
        error("Fuselage can only be closed to a 'point' or an 'edge' but '$fuse_end' was provided.")
    else
        fuselage.layout.opt_tapers_to = fuse_end
    end

    fuselage.layout.x_nose = Distance(readgeom("x_nose_tip")) 
    fuselage.layout.x_pressure_shell_fwd = Distance(readgeom("x_pressure_shell_fwd"))
    fuselage.layout.x_pressure_shell_aft = Distance(readgeom("x_pressure_shell_aft"))
    fuselage.layout.x_start_cylinder = Distance(readgeom("x_start_cylinder"))
    fuselage.layout.x_end_cylinder = Distance(readgeom("x_end_cylinder"))
    fuselage.layout.x_cone_end = Distance(readgeom("x_cone_end"))
    fuselage.layout.x_end = Distance(readgeom("x_end")) 
    fuselage.layout.l_cabin_cylinder = fuselage.layout.x_end_cylinder - fuselage.layout.x_start_cylinder


aero = read_input("Aero", fuse, dfuse)
daero = dfuse["Aero"]
readaero(x) = read_input(x, aero, daero)
    para[iafexcdf, :, :] .= transpose(readaero("excrescence_drag_factor")) #transpose for proper vector broadcasting
    para[iafduo, :, :] .= transpose(readaero("wingroot_fuse_overspeed"))
    para[iafdus, :, :] .= transpose(readaero("wingbreak_fuse_overspeed"))
    para[iafdut, :, :] .= transpose(readaero("wingtip_fuse_overspeed"))

    # if calculates_pitching_moment_volume, 
    #  CMVf1 will be calculated at sizing using slender body theory assumptions
    # else, the value specified in the input is used. 
    # Note: specified CMVf1 is only accessible when read_aircraft_model() is called.
    fuselage.calculates_pitching_moment_volume = readaero("calculates_pitching_moment_volume")
    parg[igCMVf1] = fuselage.calculates_pitching_moment_volume ? nothing : Vol(readaero("pitching_moment_volume"))
    #in either case, the CL where Mfuse = 0 must be specified
    parg[igCLMf0] = readaero("CL_zero_fuse_moment")
    
    parg[igfBLIf] = readaero("BLI_frac") #fuselage boundary layer ingestion fraction

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

    fuselage.APU.W = readweight("APU_weight_fraction")*exitlimit*Wpax
    fuselage.seat.W = readweight("seat_weight_fraction")*exitlimit*Wpax
    fuselage.added_payload.W = readweight("add_payload_weight_fraction")*exitlimit*Wpax

# ------ End fuse -------

# ---------------------------------
# Landing gear
# ---------------------------------
lg = read_input("LandingGear", data, default)
dlg = default["LandingGear"]
readlg(x::String) = read_input(x, lg, dlg)

#Landing gear CG positions or offsets
xlgnose = Distance(readlg("x_nose_landing_gear"))
landing_gear.nose_gear.weight = TASOPT.structures.Weight(W = 0.0, x = xlgnose)
landing_gear.main_gear.distance_CG_to_landing_gear = Distance(readlg("x_main_landing_gear_offset"))

#The mass model for the landing gear can be specified by the user
lgmodel = readlg("landing_gear_model")
landing_gear.model = lgmodel

if lowercase(lgmodel) == "mass_fractions" #This is the most basic model, just fixed fractions of the MTOW
    landing_gear.nose_gear.overall_mass_fraction = readlg("LG_nose_weight_fraction")
    landing_gear.main_gear.overall_mass_fraction = readlg("LG_main_weight_fraction")
elseif lowercase(lgmodel) == "historical_correlations" #model based on historical-data relations in Raymer (2012)
    landing_gear.main_gear.y_offset_halfspan_fraction = readlg("y_main_landing_gear_halfspan_fraction")
    landing_gear.tailstrike_angle = Angle(readlg("tailstrike_angle"))
    landing_gear.wing_dihedral_angle = Angle(readlg("wing_dihedral_angle")) #TODO consider storing this as a wing parameter
    landing_gear.engine_ground_clearance = Distance(readlg("engine_ground_clearance"))
    landing_gear.nose_gear.number_struts = readlg("LG_nose_number_struts")
    landing_gear.nose_gear.wheels_per_strut = readlg("LG_nose_wheels_per_strut")
    landing_gear.main_gear.number_struts = readlg("LG_main_number_struts")
    landing_gear.main_gear.wheels_per_strut = readlg("LG_main_wheels_per_strut")
end
# ------ End landing gear -------

#Fuel storage options
fuse_tank = fuselage_tank() #Initialize struct for fuselage fuel tank params

has_wing_fuel = readfuel("fuel_in_wing")
if !(has_wing_fuel) #If fuel is stored in fuselage
    fuel_stor = readfuel("Storage")
    dfuel_stor = dfuel["Storage"]
    readfuel_storage(x::String) = read_input(x, fuel_stor, dfuel_stor)

    fuse_tank.placement = readfuel_storage("tank_placement")
    fuse_tank.fueltype = fueltype
    fuse_tank.clearance_fuse = Distance(readfuel_storage("fuselage_clearance"))

    fuse_tank.sizes_insulation = readfuel_storage("sizes_insulation")
    fuse_tank.t_insul = readfuel_storage("insulation_segment_base_thickness")
    insul_mats_names = readfuel_storage("insulation_material")
    insul_mats = []
    for insul_mat_name in insul_mats_names
        push!(insul_mats, ThermalInsulator(insul_mat_name))
    end
    fuse_tank.material_insul = insul_mats
    if fuse_tank.sizes_insulation
        fuse_tank.boiloff_rate = readfuel_storage("cruise_boiloff_rate")
        fuse_tank.iinsuldes = readfuel_storage("insulation_thicknesses_design_indices")
    end
    
    inner_mat_name = readfuel_storage("inner_vessel_material")
    fuse_tank.inner_material = StructuralAlloy(inner_mat_name)
    
    fuse_tank.ARtank = readfuel_storage("tank_aspect_ratio")
    fuse_tank.theta_inner = Angle(readfuel_storage("inner_vessel_support_angle"))

    fuse_tank.pvent = Pressure(readfuel_storage("pressure_venting"))
    fuse_tank.pinitial = Pressure(readfuel_storage("pressure_initial"))
    fuse_tank.t_hold_orig = Times(readfuel_storage("hold_departure"))
    fuse_tank.t_hold_dest = Times(readfuel_storage("hold_arrival"))
    
    fuse_tank.ftankadd = readfuel_storage("additional_mass_fraction")
    fuse_tank.ew = readfuel_storage("weld_efficiency")
    fuse_tank.ullage_frac = readfuel_storage("ullage_fraction")
    fuse_tank.qfac = readfuel_storage("heat_leak_factor")
    fuse_tank.pfac = readfuel_storage("pressure_rise_factor")

    #Store takeoff temperatures in tank object as well for ease of access
    for (i,altTO) in enumerate(parm[imaltTO, :]/1e3)
        T_std, _, _, _, _ = atmos(altTO)
        push!(fuse_tank.TSLtank, parm[imT0TO,i] - T_std + Tref)
    end
    
    has_vacuum = TASOPT.CryoTank.check_vacuum(fuse_tank.material_insul) #flag to check if an outer vessel is needed

    if has_vacuum #If tank is double-walled
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
        nftanks = 1
    elseif (fuse_tank.placement == "both") 
        nftanks = 2
    end
    fuse_tank.tank_count = nftanks
else #else all fuel in wings
    nftanks = 0
end #if

# ---------------------------------
# Wing
# ---------------------------------
# Setup wing
wing_i = read_input("Wing", data, default)
dwing = default["Wing"]
readwing(x) = read_input(x, wing_i, dwing)
    wing.has_strut = readwing("has_strut")

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

    ## Strut details only used if has_strut is true
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
    htail.tip_lift_loss = wing.tip_lift_loss
    vtail.tip_lift_loss = wing.tip_lift_loss

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
    #TODO: change 
    multi_section = readhtail("multi_section")
    # if !multi_section
    #     htail.layout.ηs = 0.0
    # end
    htail.inboard.λ = 1.0
    
    # igbs = igbo
    # strutz = 0
    # lambdat = gammat = iglambdah 
    # lambdas = gammas = 1.0

    # create inner
    # lambdas, gammas = 1.0
    # igbs = igbo
    # hboxs = hboxh
    htail.layout.sweep = readhtail("sweep")
    htail.layout.root_span = 2*Distance(readhtail("center_box_halfspan"))

    htail.layout.box_x  = Distance(readhtail("x_Htail"))
    htail.layout.z = Distance(readhtail("z_Htail"))
    htail.ntails  = readhtail("number_Htails")

    htail.CL_CLmax = readhtail("max_tail_download")

    htail_sizing = readhtail("opt_sizing")
    if compare_strings(htail_sizing,"fixed_Vh")
        htail.opt_sizing = htail_sizing
        htail.volume = readhtail("Vh")
    elseif compare_strings(htail_sizing,"CLmax_fwdCG")
        htail.opt_sizing = htail_sizing
        htail.CL_max_fwd_CG = readhtail("CLh_at_max_forward_CG")
        htail.volume = 1.0
    else
        error("Horizontal tail can only be sized via:
            \"fixed_Vh\":   specified tail volume coeff (\"Vh\");
            \"CLmax_fwdCG\": specified CLh (\"CLh_at_max_forward_CG\") at 'worst-case': max-forward CG, max wing lift")
    end

    opt_move_wing = readhtail("opt_move_wing")
    #if not an expected input, throw an error
    valid_options = ["fixed", "fixed_CLh", "min_static_margin"]
    if !(any(compare_strings.(opt_move_wing, valid_options)))
        error("Input error: \"opt_move_wing\" = $opt_move_wing\nWing position during horizontal tail sizing can only be sized via:\n" *
              join([">\"$opt\"" for opt in valid_options], "\n"))
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

    htail.inboard.cross_section.width_to_chord = readhtail("box_width_to_chord")
    htail.inboard.cross_section.thickness_to_chord = readhtail("box_height_chord")
    htail.inboard.cross_section.web_to_box_height  = readhtail("web_height_hbox")
    

vtail_input = readtails("Vtail")
dvtail = dtails["Vtail"]
readvtail(x) = read_input(x, vtail_input, dvtail)
    vtail.layout.AR = readvtail("AR_Vtail")
    vtail.outboard.λ = readvtail("taper")
    vtail.inboard.λ = 1.0
    vtail.layout.sweep  = readvtail("sweep")
    vtail.layout.root_span = Distance(readvtail("center_box_halfspan"))
    vtail.layout.box_x  = Distance(readvtail("x_Vtail"))
    vtail.ntails  = readvtail("number_Vtails")

    vtail_sizing = readvtail("opt_sizing")
    if compare_strings(vtail_sizing, "fixed_Vv")
        vtail.opt_sizing = vtail_sizing
        vtail.volume = readvtail("Vv")
    elseif compare_strings(vtail_sizing, "OEI")
        vtail.opt_sizing = vtail_sizing
        parg[igCLveout] = readvtail("CLv_at_engine_out")
    else
        error("Vertical tail can only be sized via:
            \"fixed_Vv\": specified tail volume coeff \"Vv\";
            \"OEI\": specified CL at one engine out trim (\"OEI\") via \"CLv_at_engine_out\"")
    end

    vtail.CL_max = readvtail("CLv_max")

    vtail.weight_fraction_added = readvtail("added_weight_fraction")
    vtail.outboard.cross_section.width_to_chord = readvtail("box_width_to_chord")
    vtail.outboard.cross_section.thickness_to_chord = readvtail("box_height_chord")
    vtail.outboard.cross_section.web_to_box_height  = readvtail("web_height_hbox")

    vtail.inboard.cross_section.width_to_chord = readvtail("box_width_to_chord")
    vtail.inboard.cross_section.thickness_to_chord = readvtail("box_height_chord")
    vtail.inboard.cross_section.web_to_box_height  = readvtail("web_height_hbox")

# ----- End Stabilizers -----

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
propsys = read_input("prop_sys_arch", options, doptions)
prop = read_input("Propulsion", data, default)
dprop = default["Propulsion"]
readprop(x) = read_input(x, prop, dprop)

parg[igneng] = readprop("number_of_engines")

if lowercase(propsys) == "tf"
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
    eng_has_BLI_cores = !readprop("core_in_clean_flow")

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

    #HPT cooled efficiency
    pare[iedehtdfc,:,:] .= readcool("HPT_efficiency_derivative_with_cooling")
    pare[iefc0,:,:] .= readcool("baseline_cooling_fraction")

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

elseif lowercase(propsys) == "constant_tsfc" #For constant TSFC model
    ROCdes = readprop("rate_of_climb")
    if ROCdes isa AbstractVector
        para[iaROCdes,ipclimb1:ipclimbn,:] .= [Speed(x) for x in ROCdes]
    else
        para[iaROCdes,ipclimb1:ipclimbn,:] .= Speed(ROCdes)
    end
    pare[ieTSFC,ipclimb1:ipclimbn,:] .= readprop("climb_TSFC")
    pare[ieTSFC,ipcruise1:ipcruisen,:] .= readprop("cruise_TSFC")
    pare[ieTSFC,ipdescent1:ipdescentn,:] .= readprop("descent_TSFC")

else #unrecognized input
    @warn("The engine type is not recognized")
end

nac = readprop("Nacelles")
dnac = dprop["Nacelles"]
    #- nacelle drag stuff
    parg[igrSnace] = read_input("nacelle_pylon_wetted_area_ratio", nac, dnac)
    parg[igrVnace] = read_input("nacelle_local_velocity_ratio", nac, dnac)

weight = readprop("Weight")
dweight = dprop["Weight"]
    parg[igfeadd] = read_input("engine_access_weight_fraction", weight, dweight)
    parg[igfpylon] = read_input("pylon_weight_fraction", weight, dweight)
    
    #read/check engine weight model options
    TF_wmodel = read_input("weight_model", weight, dweight)
    if compare_strings(propsys, "tf")
    #TODO: reincorporate "pantalone_basic" and "pantalone_adv" for direct-drive turbofans
        engineweightname = TF_wmodel
        engineweight! = tfweightwrap!
        if !(TF_wmodel in ["md", "fitzgerald_basic", "fitzgerald_adv"]) 
            error("\"$TF_wmodel\" engine weight model was specifed. 
            Engine weight can only be \"MD\", \"fitzgerald_basic\" or \"fitzgerald_adv\".")
        end
    elseif compare_strings(TF_wmodel, "fractional_weight")
        parg[igfeng] = read_input("engine_weight_fraction", weight, dweight)
        engineweightname = "fractional_weight"
        engineweight! = TASOPT.engine.fractional_engine_weight!
    elseif compare_strings(TF_wmodel, "constant_weight")
        eng_weight = Force(read_input("engine_weight_total", weight, dweight))
        engineweightname = "constant_weight"
        function constant_engine_weight(ac)
            ac.parg[igWeng] = eng_weight
        end
        engineweight! = constant_engine_weight
    elseif compare_strings(propsys, "te")
        @warn("Propulsion weight models for turboelectric are currently not available.")
    end

# Create engine object
if compare_strings(propsys,"tf")
    modelname = "turbofan_md"
    enginecalc! = tfwrap!


    enginemodel = TASOPT.engine.TurbofanModel(modelname, enginecalc!, engineweightname, engineweight!, eng_has_BLI_cores)
    engdata = TASOPT.engine.EmptyData()
elseif compare_strings(propsys,"fuel_cell_with_ducted_fan")
    modelname = lowercase(propsys)
    engineweightname = "nasa"

    enginecalc! = calculate_fuel_cell_with_ducted_fan!
    engineweight! = fuel_cell_with_ducted_fan_weight!
    enginemodel = TASOPT.engine.FuelCellDuctedFan(modelname, enginecalc!, engineweightname, engineweight!, eng_has_BLI_cores)
    pare[iePfanmax,:,:] .= 20e6

    fcdata = TASOPT.engine.FuelCellDuctedFanData(2)

    fcdata.type = "HT-PEMFC"
    fcdata.current_density[iprotate,:] .= 1e4
    fcdata.FC_temperature .= 453.15
    fcdata.FC_pressure .= 3e5
    fcdata.water_concentration_anode .= 0.1
    fcdata.water_concentration_cathode .= 0.1
    fcdata.λ_H2 .= 3.0
    fcdata.λ_O2 .= 3.0
    fcdata.thickness_membrane = 100e-6
    fcdata.thickness_anode  = 250e-6
    fcdata.thickness_cathode  = 250e-6
    fcdata.design_voltage = 200.0
    pare[ieRadiatorepsilon,:,:] .= 0.7
    pare[ieRadiatorMp,:,:] .= 0.12
    pare[ieDi,:,:] .= 0.4

    para[iaROCdes, ipclimb1:ipclimbn,:] .= 500 * ft_to_m / 60
    engdata = fcdata

elseif lowercase(propsys) == "constant_tsfc"
    modelname = "constant_TSFC"
    enginecalc! = TASOPT.engine.constant_TSFC_engine!
    
    enginemodel = TASOPT.engine.TurbofanModel(modelname, enginecalc!, engineweightname, engineweight!, false)
    calculate_takeoff = false #Engine model cannot be used for takeoff

else
    #TODO: FIX THIS labelling
    error("Propulsion system performance model \"$propsys\" specified. Choose between
    > `TF` - detailed turbofan cycle model
    > `constant_tsfc` - simplest engine model
    > `fuel_cell_with_ducted_fan` - hydrogen fuel cell driving electric fans")
end
engine = TASOPT.engine.Engine(enginemodel, engdata, Vector{TASOPT.engine.HeatExchanger}())
    
# Heat exchangers, only if in the model file (assigned 0 in the default .tomls)
# Would be better if there were an independent toggle like `has_wing_fuel` for `fuse_tank` but good enough for now
if "HeatExchangers" in keys(prop) && !isempty(prop["HeatExchangers"])
    HEx = readprop("HeatExchangers")
    dHEx = dprop["HeatExchangers"]
    HX_add_mass = read_input("added_mass_frac", HEx, dHEx)
    
    has_recirculation = read_input("has_recirculation", HEx, dHEx)
    recircT = Temp(read_input("recirculation_temperature", HEx, dHEx))
    pare[ieDi, :, :] .= Distance(read_input("core_inner_diameter", HEx, dHEx))
    HXmaxL = Distance(read_input("maximum_heat_exchanger_length", HEx, dHEx))
    
    PreCorder = read_input("precooler_order", HEx, dHEx)
    PreCepsilon = read_input("precooler_effectiveness", HEx, dHEx)
    PreCMp = read_input("precooler_inlet_mach", HEx, dHEx)

    InterCorder = read_input("intercooler_order", HEx, dHEx)
    InterCepsilon = read_input("intercooler_effectiveness", HEx, dHEx)
    InterCMp = read_input("intercooler_inlet_mach", HEx, dHEx)

    Regenorder = read_input("regenerative_order", HEx, dHEx)
    Regenepsilon = read_input("regenerative_effectiveness", HEx, dHEx)
    RegenMp = read_input("regenerative_inlet_mach", HEx, dHEx)

    TurbCorder = read_input("turbine_cooler_order", HEx, dHEx)
    TurbCepsilon = read_input("turbine_cooler_effectiveness", HEx, dHEx)
    TurbCMp = read_input("turbine_cooler_inlet_mach", HEx, dHEx)

    all_eps = [PreCepsilon, InterCepsilon, Regenepsilon, TurbCepsilon]
    all_types = ["PreC", "InterC", "Regen", "TurbC"]
    all_orders = [PreCorder, InterCorder, Regenorder, TurbCorder]
    all_Mp = [PreCMp, InterCMp, RegenMp, TurbCMp]

    HXtypes = []
    Mp_in = []
    ε_des = []
    sort_i = sortperm(all_orders) #Sort according to order

    for ind in sort_i
        if (all_eps[ind] > 0) && (all_eps[ind] <= 1) #If effectiveness is between 0 and 1
                push!(HXtypes, all_types[ind])
                push!(Mp_in, all_Mp[ind])
                push!(ε_des, all_eps[ind])
        end
    end
    nmis = size(pare)[3]
            
    for (i,HXtype) in enumerate(HXtypes)
        #Create heat exchanger struct
        HX = TASOPT.engine.make_HeatExchanger(nmis)
        HX.type = HXtype
        HX.order = i
        HX.design_effectiveness = ε_des[i]
        HX.design_Mach = Mp_in[i]
        HX.added_mass_fraction = HX_add_mass
        HX.maximum_length = HXmaxL
        if i == 1 #For first HX
            HX.has_recirculation = has_recirculation
            HX.recirculation_temperature = recircT
        end
        #Add to engine heat exchangers
        push!(engine.heat_exchangers, HX)
    end
end

#create options object
ac_options = TASOPT.Options(
    opt_fuel = fueltype,
    ifuel = ifuel,
    has_wing_fuel = has_wing_fuel,
    has_centerbox_fuel = has_centerbox_fuel,
    has_fuselage_fuel = (nftanks>0),
    
    opt_engine_location = engloc,
    opt_prop_sys_arch = propsys,
    calculate_takeoff = calculate_takeoff,
    
    is_doubledecker = is_doubledecker,

    opt_move_wing = opt_move_wing,

    trefftz_config = TREFFTZ_CONFIG
)
    
#Create aircraft object
ac = TASOPT.aircraft(name, description, ac_options,
    parg, parm, para, pare, is_sized, 
    fuselage, fuse_tank, wing, htail, vtail, engine, landing_gear)

# ---------------------------------
# Recalculate cabin length
if calculate_cabin #Resize the cabin if desired, keeping deltas
    @info "Fuselage and stabilizer layouts have been overwritten; deltas will be maintained."
    update_fuse_for_pax!(ac) #update fuselage dimensions
end

return ac

end

function load_default_model()
    println("Loading default aircraft model")
    read_aircraft_model()
end