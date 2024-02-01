using TOML

"""
    read_input(k::String, dict::AbstractDict=data, 
    default_dict::AbstractDict = default)

Reads the input from a given dictonary (typically parsed from a TOML file).
If requested input does not exist in dictonary, looks for value in default input
and stores default value into the given dictonary (primarily for later output/
saving as an `aircraft` model file)
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
            input file to see all available input options.")
        end
    end
end

# Convinence functions to convert to SI units
Speed(x)    = convertSpeed(parse_unit(x)...)
Len(x)      = convertDist(parse_unit(x)...)
Force(x)    = convertForce(parse_unit(x)...)
Pressure(x) = convertPressure(parse_unit(x)...)
Stress = Pressure
Density(x)  = convertDensity(parse_unit(x)...)
Area(x)     = convertArea(parse_unit(x)...)
Vol(x)      = convertVolume(parse_unit(x)...)
Angle(x)    = convertAngle(parse_unit(x)...)
Temp(x)     = convertTemp(parse_unit(x)...)

"""
    read_aircraft_model(datafile; 
    defaultfile = joinpath(TASOPT.__TASOPTroot__, "IO/default_input.toml"))

Reads a specified TOML file that describes a TASOPT `aircraft` model 
with a fall back to the default `aircraft` definition 
provided in \"src/IO/default_input.toml\""

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
    1: TF - turbo-fan
    2: TE - turbo-electric" )
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
        1: Engines on \"wing\"
        2: Engines on \"fuselage\"")
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
else
    error("Check fuel type")
end
pari[iifwcen]  = readfuel("fuel_in_wingcen")
parg[igrWfmax] = readfuel("fuel_usability_factor")
pare[ieTft, :, :] .= readfuel("fuel_temp") #Temperature of fuel in fuel tank
pare[ieTfuel, :, :] .= readfuel("fuel_temp") #Initialize fuel temperature as temperature in tank
parg[igrhofuel] = readfuel("fuel_density")

#Fuel storage options
fuse_tank = fuselage_tank() #Initialize struct for fuelage fuel tank params

fuel_stor = readfuel("Storage")
dfuel_stor = dfuel["Storage"]
readfuel_storage(x::String) = read_input(x, fuel_stor, dfuel_stor)

pari[iifwing]  = readfuel_storage("fuel_in_wing")

if pari[iifwing]  == 0 #If fuel is stored in fuselage
    pari[iinftanks] = readfuel_storage("fuel_tanks_in_fuse")

    fuse_tank.t_insul = readfuel_storage("insulation_segment_base_thickness")
    fuse_tank.k_insul= readfuel_storage("insulation_segment_conductivity")
    fuse_tank.rho_insul = readfuel_storage("insulation_segment_density")
    fuse_tank.iinsuldes = readfuel_storage("insulation_thicknesses_design_indices")
    fuse_tank.sigskin = readfuel_storage("skin_yield_strength")
    fuse_tank.rhoskintank = readfuel_storage("tank_skin_density")
    fuse_tank.max_boiloff = readfuel_storage("maximum_boiloff_rate")
    fuse_tank.ARtank = readfuel_storage("fuel_tanks_in_fuse")
    fuse_tank.clearance_fuse = readfuel_storage("fuselage_clearance")
    fuse_tank.ptank = readfuel_storage("tank_pressure")
    fuse_tank.ftankstiff = readfuel_storage("stiffener_mass_fraction")
    fuse_tank.ftankadd = readfuel_storage("additional_mass_fraction")
end

# Setup mission variables
ranges = readmis("range")
parm[imRange, :] .= Len.(ranges)

Wpax =  Force(readmis("weight_per_pax"))
parm[imWpay, :] .= readmis("pax") * Wpax
parg[igWpaymax] = readmis("max_pax") * Wpax
parg[igfreserve] = readmis("fuel_reserves")
parg[igVne] = Speed(readmis("Vne"))
parg[igNlift] = readmis("Nlift")

##Takeoff
takeoff = readmis("Takeoff")
dtakeoff = dmis["Takeoff"]
readtakeoff(x) = read_input(x, takeoff, dtakeoff)
parm[imaltTO, :] .= Len.(readtakeoff("takeoff_alt"))
parg[igmubrake] = readtakeoff("braking_resistance_coeff")
parg[igmuroll]  = readtakeoff("rolling_resistance_coeff")
parg[ighobst]   = Len(readtakeoff("takeoff_obstacle_height"))
parg[igcdefan]  = readtakeoff("CD_dead_engine")
parg[igCDgear]  = readtakeoff("CD_landing_gear")
parg[igCDspoil] = readtakeoff("CD_spoilers")
parg[iglBFmax]  = Len(readtakeoff("max_balanced_field_length"))
parg[igNland]   = readtakeoff("Nland")

T0TO = Temp.(readtakeoff("takeoff_T"))
parm[imT0TO, :] .= T0TO 
para[iaclpmax, ipclimb1, :] .= readtakeoff("CL_max_perp")
para[iaclpmax, ipstatic:ipcutback, :] .= readtakeoff("CL_max_perp")
para[iaclpmax, ipdescentn, :] .= readtakeoff("CL_max_perp")

##Climb parameters
climb = readmis("Climb")
dclimb = dmis["Climb"]
parg[iggtocmin] = Angle(read_input("minimum_top-of-climb_gradient",
                 climb, dclimb))

##Cruise parameters
cruise = readmis("Cruise")
dcruise = dmis["Cruise"]
readcruise(x) = read_input(x, cruise, dcruise)
para[iaalt, ipcruise1, :] .= Len.(readcruise("cruise_alt"))
para[iaMach, ipclimbn:ipdescent1, :] .= readcruise("cruise_mach")
para[iaCL, ipclimb1+1:ipdescentn-1, :] .= readcruise("cruise_CL")

##Descent parameters
des = readmis("Descent")
ddes = dmis["Descent"]
readdes(x) = read_input(x, des, ddes)
parm[imgamVDE1, :] .= Angle(readdes("descent_angle_top-of-descent"))
parm[imgamVDEn, :] .= Angle(readdes("descent_angle_bottom-of-descent"))

#---------- End Mission vars --------------

# ---------------------------------
# Fuselage
# ---------------------------------
# Setup Fuselage 
fuse = read_input("Fuselage", data, default)
dfuse = default["Fuselage"]
cabinPressureAlt_km = convertDist(parse_unit(read_input("cabin_pressure_altitude",
                                            fuse, dfuse))..., "km")
_, p_cabin, _, _, _ = atmos(cabinPressureAlt_km)
parg[igpcabin] = p_cabin

aero = read_input("Aero", fuse, dfuse)
daero = dfuse["Aero"]
readaero(x) = read_input(x, aero, daero)
    para[iafexcdf, 1:iptotal, :] .= readaero("excrescence_drag_factor")

    para[iafduo, :, :] .= readaero("wingroot_fuse_overspeed")
    para[iafdus, :, :] .= readaero("wingbreak_fuse_overspeed")
    para[iafdut, :, :] .= readaero("wingtip_fuse_overspeed")

    parg[igCMVf1] = Vol(readaero("fuse_moment_volume_deriv"))
    parg[igCLMf0] = readaero("CL_zero_fuse_moment")
    
    parg[igfBLIf] = readaero("BLI_frac")

weight = read_input("Weights", fuse, dfuse)
dweight = dfuse["Weights"]
readweight(x) = read_input(x, weight, dweight)
    parg[igfframe]  = readweight("frame")
    parg[igfstring] = readweight("stringer")
    parg[igffadd]   = readweight("additional")

    parg[igWfix] = Force(readweight("fixed_weight"))

    parg[igWpwindow] = readweight("window_per_length")
    parg[igWppinsul] = readweight("window_insul_per_area")
    parg[igWppfloor] = readweight("floor_weight_per_area")

    parg[igfhpesys] = readweight("HPE_sys_weight_fraction")
    parg[igflgnose] = readweight("LG_nose_weight_fraction")
    parg[igflgmain] = readweight("LG_main_weight_fraction")

    parg[igfapu] = readweight("APU_weight_fraction")
    parg[igfseat] = readweight("seat_weight_fraction")
    parg[igfpadd] = readweight("add_payload_weight_fraction")

geom = read_input("Geometry", fuse, dfuse)
dgeom = dfuse["Geometry"]
readgeom(x) = read_input(x, geom, dgeom)
    parg[igRfuse]  = Len(readgeom("radius"))
    parg[igdRfuse] = Len(readgeom("dRadius"))
    parg[igwfb]    = Len(readgeom("y_offset"))
    parg[ighfloor] = Len(readgeom("floor_depth"))
    parg[ignfweb]  = readgeom("Nwebs")

    parg[iganose] = readgeom("a_nose")
    parg[igbtail] = readgeom("b_tail")

    fuse_end = lowercase(readgeom("taper_fuse_to")) 
    if fuse_end == "point"
        pari[iifclose] = 0
    elseif fuse_end == "edge"
        pari[iifclose] = 1
    else
        pari[iifclose] = 0
        @warn "Fuselage can only be closed to a 'point' or an 'edge'"*
                " but '$fuse_end' was provided."*
                " Setting fuselage to end at a point."
    end

    parg[iglambdac] = readgeom("tailcone_taper")
    parg[igrMh] = readgeom("HT_load_fuse_bend_relief")
    parg[igrMv] = readgeom("VT_load_fuse_bend_relief")

    parg[igxnose]   = Len(readgeom("x_nose_tip")) 
    parg[igxshell1] = Len(readgeom("x_pressure_shell_fwd"))
    parg[igxblend1] = Len(readgeom("x_start_cylinder"))
    parg[igxblend2] = Len(readgeom("x_end_cylinder"))
    parg[igxshell2] = Len(readgeom("x_pressure_shell_aft"))
    parg[igxconend] = Len(readgeom("x_cone_end"))
    parg[igxend]    = Len(readgeom("x_end")) 
    
    parg[igxlgnose]  = Len(readgeom("x_nose_landing_gear"))
    parg[igdxlgmain] = Len(readgeom("x_main_landing_gear_offset"))
    parg[igxapu]     = Len(readgeom("x_APU"))
    parg[igxhpesys]  = Len(readgeom("x_HPE_sys"))

    parg[igxfix] = Len(readgeom("x_fixed_weight"))

    parg[igxeng] = Len(readgeom("x_engines"))
    parg[igyeng] = Len(readgeom("y_critical_engines"))

    parg[iglblend2blend] = parg[igxblend2] - parg[igxblend1]

# ------ End fuse -------
# ---------------------------------
# Wing
# ---------------------------------
# Setup wing
wing = read_input("Wing", data, default)
dwing = default["Wing"]
readwing(x) = read_input(x, wing, dwing)
    pari[iiwplan] = readwing("wing_planform")
    if readwing("strut_braced_wing")
        pari[iiwplan] = 2
    end

    parg[igsweep] = readwing("sweep")
    parg[igAR] = readwing("AR")
    parg[igbmax] = Len(readwing("maxSpan"))

    parg[iglambdas] = readwing("inner_panel_taper_ratio")
    parg[iglambdat] = readwing("outer_panel_taper_ratio")
    parg[igetas]    = readwing("panel_break_location")

    parg[igbo] = 2*Len(readwing("center_box_halfspan"))
    parg[igwbox]  = readwing("box_width_chord")
    parg[ighboxo] = readwing("root_thickness_to_chord")
    parg[ighboxs] = readwing("spanbreak_thickness_to_chord")
    parg[igrh]    = readwing("hweb_to_hbox")
    parg[igXaxis] = readwing("spar_box_x_c")

    parg[igxwbox] = Len(readwing("x_wing_box"))
    parg[igzwing] = Len(readwing("z_wing"))

    ## Strut details only used if strut_braced_wing is true
    parg[igzs]      = Len(readwing("z_strut"))
    parg[ighstrut]  = readwing("strut_toc")
    parg[igrVstrut] = readwing("strut_local_velocity_ratio")

# ----------------------------------
# ------- Wing Aerodynamics --------
# ----------------------------------
aero = readwing("Aero")
daero = dwing["Aero"]
    parg[igfLo] = readaero("fuselage_lift_carryover_loss_factor")
    parg[igfLt] = readaero("wing_tip_lift_rolloff_factor")

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
    para[iarcls, 1:ipclimb1, :] .= readtakeoff("cls_clo")    #  rcls    break/root cl ratio = cls/clo
    para[iarclt, 1:ipclimb1, :] .= readtakeoff("clt_clo")    #  rclt    tip  /root cl ratio = clt/clo
    para[iacmpo, 1:ipclimb1, :] .= readtakeoff("cm_o")      #  cmpo    root  cm
    para[iacmps, 1:ipclimb1, :] .= readtakeoff("cm_s")      #  cmps    break cm
    para[iacmpt, 1:ipclimb1, :] .= readtakeoff("cm_t")      #  cmpt    tip   cm

# Clean climb cruise descent and for wing structure sizing
climb = readaero("Climb")
dclimb = daero["Climb"]
readclimb(x) = read_input(x, climb, dclimb)
    para[iarcls, ipclimb1+1:ipdescentn-1, :] .= readclimb("cls_clo")   #  rcls    break/root cl ratio = cls/clo
    para[iarclt, ipclimb1+1:ipdescentn-1, :] .= readclimb("clt_clo")   #  rclt    tip  /root cl ratio = clt/clo
    para[iacmpo, ipclimb1+1:ipdescentn-1, :] .= readclimb("cm_o")      #  cmpo    root  cm
    para[iacmps, ipclimb1+1:ipdescentn-1, :] .= readclimb("cm_s")      #  cmps    break cm
    para[iacmpt, ipclimb1+1:ipdescentn-1, :] .= readclimb("cm_t")      #  cmpt    tip   cm

# Landing, forward CG tail sizing case
land = readaero("Landing")
dland = daero["Landing"]
readland(x) = read_input(x, land, dland)
    para[iarcls, ipdescentn, :] .= readland("cls_clo")   #  rcls    break/root cl ratio = cls/clo
    para[iarclt, ipdescentn, :] .= readland("clt_clo")   #  rclt    tip  /root cl ratio = clt/clo
    para[iacmpo, ipdescentn, :] .= readland("cm_o")      #  cmpo    root  cm
    para[iacmps, ipdescentn, :] .= readland("cm_s")      #  cmps    break cm
    para[iacmpt, ipdescentn, :] .= readland("cm_t")      #  cmpt    tip   cm
    
# ----------------------------------
# ---------- Wing Weight -----------
# ----------------------------------
# Wing weight fractions of flight surfaces
# and secondary wing components, 
weight = readwing("Weightfracs")
dweight = dwing["Weightfracs"]

    parg[igfflap] = readweight("flap")
    parg[igfslat] = readweight("slat")
    parg[igfaile] = readweight("aileron")
    parg[igflete] = readweight("leading_trailing_edge")
    parg[igfribs] = readweight("ribs")
    parg[igfspoi] = readweight("spoilers")
    parg[igfwatt] = readweight("attachments")

# ---- End Wing -----

# ---------------------------------
# Stabilizers
# ---------------------------------
tails = read_input("Stabilizers", data, default)
dtails = default["Stabilizers"]
readtails(x) = read_input(x, tails, dtails)
    para[iacdft, 1:iptotal, :]   .= readtails("lowspeed_cdf")  #  cdft    tail profile cd
    para[iacdpt, 1:iptotal, :]   .= readtails("lowspeed_cdp")  #  cdpt    
    para[iaRereft, 1:iptotal, :] .= readtails("Re_ref")  #  Rereft  

    para[iafexcdt, 1:iptotal, :] .= readtails("excrescence_drag_factor")

    htail = readtails("Htail")
    dhtail = dtails["Htail"]

readhtail(x) = read_input(x, htail, dhtail)
    parg[igARh]     = readhtail("AR_Htail")
    parg[iglambdah] = readhtail("taper")
    parg[igsweeph]  = readhtail("sweep")
    parg[igboh]     = 2*Len(readhtail("center_box_halfspan"))

    parg[igxhbox]  = Len(readhtail("x_Htail"))
    parg[igzhtail] = Len(readhtail("z_Htail"))

    parg[igCLhNrat] = readhtail("max_tail_download")

    htail_size = lowercase(readhtail("HTsize"))
    if htail_size == "vh"
        pari[iiHTsize] = 1
        parg[igVh] = readhtail("Vh")
    elseif htail_size == "maxforwardcg"
        pari[iiHTsize] = 2
        parg[igCLhCGfwd] = readhtail("CLh_at_max_forward_CG")
        parg[igVh] = 1.0
    else
        error("Horizontal tail can only be sized via:
            1: specified tail volume coeff \"Vh\";
            2: specified CLh at max-forward CG case during landing (\"maxforwardCG\")")
    end


    movewing = readhtail("move_wingbox")
    if typeof(movewing) == Int
        pari[iixwmove] = movewing
    elseif typeof(movewing) <: AbstractString
        movewing = lowercase(movewing)
        if movewing =="fix"
            pari[iixwmove] = 0
        elseif movewing == "clhspec"
            pari[iixwmove] = 1
        elseif movewing == "smmin"
            pari[iixwmove] = 2
        else
            error("Wing position duirng horizontal tail sizing can only be sized via:
            0: \"fix\" wing position;
            1: move wing to get CLh=\"CLhspec\" in cruise 
            2: move wing to get min static margin = \"SMmin\"")
        end
    else
        error("Check wing position input during htail sizing... something isn't right")
    end
    parg[igSMmin] = readhtail("SM_min")

    parg[igCLhspec] = readhtail("CLh_spec")

    parg[igdepsda] = readhtail("downwash_factor")
    parg[igdCLnda] = readhtail("nacelle_lift_curve_slope")

    parg[igfCDhcen] = readhtail("CD_Htail_from_center")
    parg[igCLhmax]  = readhtail("CLh_max")

    parg[igfhadd] = readhtail("added_weight_fraction")

    parg[igwboxh] = readhtail("box_width_chord")
    parg[ighboxh] = readhtail("box_height_chord")
    parg[igrhh]   = readhtail("web_height_hbox")


vtail = readtails("Vtail")
dvtail = dtails["Vtail"]
readvtail(x) = read_input(x, vtail, dvtail)
    parg[igARv]     = readvtail("AR_Vtail")
    parg[iglambdav] = readvtail("taper")
    parg[igsweepv]  = readvtail("sweep")
    parg[igbov]     = Len(readvtail("center_box_halfspan"))
    parg[igxvbox]  = Len(readvtail("x_Vtail"))
    parg[ignvtail]  = readvtail("number_Vtails")

    vtail_size = lowercase(readvtail("VTsize"))
    if vtail_size == "vv"
        pari[iiVTsize] = 1
        parg[igVv] = readvtail("Vv")
    elseif vtail_size == "oei"
        pari[iiVTsize] = 2
        parg[igCLveout] = readvtail("CLv_at_engine_out")
    else
        error("Vertical tail can only be sized via:
            1: specified tail volume coeff \"Vv\";
            2: specified CL at one engine out trim (\"OEI\")")
    end

    parg[igCLvmax] = readvtail("CLv_max")

    parg[igfvadd] = readvtail("added_weight_fraction")
    parg[igwboxv] = readvtail("box_width_chord")
    parg[ighboxv] = readvtail("box_height_chord")
    parg[igrhv]   = readvtail("web_height_hbox")

# ----- End Stabilizers -----

# ---------------------------------
# Structures
# ---------------------------------

structures = read_input("Structures", data, default)
dstructures = default["Structures"]
readstruct(x) = read_input(x, structures, dstructures)
    parg[igsigfac] = readstruct("stress_factor")
    #- allowable stresses at sizing cases
    parg[igsigskin] = Stress(readstruct("sigma_fuse_skin"))
    parg[igsigbend] = Stress(readstruct("sigma_fuse_bending"))

    parg[igsigcap] = Stress(readstruct("sigma_caps"))
    parg[igtauweb] = Stress(readstruct("tau_webs"))

    parg[igsigstrut] = Stress(readstruct("sigma_struts"))

    #- fuselage shell modulus ratio, for bending material sizing
    parg[igrEshell] = readstruct("fuse_shell_modulus_ratio")

    #- moduli, for strut-induced buckling load estimation
    parg[igEcap] = Stress(readstruct("E_wing_spar_cap"))
    parg[igEstrut] = Stress(readstruct("E_struts"))

    #- structural material densities
    parg[igrhoskin] = Density(readstruct("skin_density"))  #  rhoskin     fuselage skin
    parg[igrhobend] = Density(readstruct("fuse_stringer_density"))  #  rhobend     fuselage bending stringers 
    parg[igrhocap]  = Density(readstruct("wing_tail_cap_density"))  #  rhocap  	wing, tail bending caps	 
    parg[igrhoweb]  = Density(readstruct("wing_tail_web_density"))  #  rhoweb  	wing, tail shear webs	 
    parg[igrhostrut]= Density(readstruct("strut_density"))  #  rhostrut	strut   

# ---------------------------------
# Propulsion systems
# ---------------------------------

prop = read_input("Propulsion", data, default)
dprop = default["Propulsion"]
readprop(x) = read_input(x, prop, dprop)
    parg[igneng] = readprop("number_of_engines")
    parg[igTmetal] = Temp(readprop("T_max_metal"))
    parg[igfTt4CL1] = readprop("Tt4_frac_bottom_of_climb")
    parg[igfTt4CLn] = readprop("Tt4_frac_top_of_climb")

    pare[ieTt4, :, :] .= Temp(readprop("Tt4_cruise"))

    Tt4TO = Temp(readprop("Tt4_takeoff"))
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
    pihc = OPR/pilc

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

    parg[igmofWpay] = mofftpax / Wpax
    parg[igmofWMTO] = mofftmMTO / gee
    parg[igPofWpay] = Pofftpax / Wpax
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

try #If heat exchanger field exists in the input file
    HEx = readprop("HeatExchangers")
    dHEx = dprop["HeatExchangers"]
        pare[iefrecirc, :, :] .= read_input("recirculation_flag", HEx, dHEx)
        pare[ierecircT, :, :] .= read_input("recirculation_temperature", HEx, dHEx)
        pare[iehlat, :, :] .= read_input("latent_heat", HEx, dHEx)
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
catch #Do nothing if the heat exchanger field does not exist
end

return TASOPT.aircraft(name, description,
pari, parg, parm, para, pare, fuse_tank)

end

function load_default_model()
    println("Loading default aircraft model")
    read_aircraft_model()
end
