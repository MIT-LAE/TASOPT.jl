"""
A trivial wrapper that runs a NPSS batch file `bat_file` at the specifed `dir`

>Note: 
>Remember to use `/` as separators
"""
function NPSS_run(dir, bat_file)
    global time_run_NPSS += @elapsed run(`cmd /c cd $dir '&&' $bat_file `)
    return nothing
end

"""
This function starts up and returns an NPSS process that can then be written to
"""
function startNPSS(dir, bat_file)
    NPSS = open(`cmd /c cd $dir '&&' $bat_file `, "w+")
    return NPSS
end

"""
Ends NPSS process that can then be written to
"""
function endNPSS(NPSS)
    write(NPSS, "999 \\n")
    if(NPSS.exitcode == 0)
        close(NPSS)
    end
end

"""
# Design mode

Writes an input file for NPSS Turboshaft model 

    - Abmient altitude and mach number
    - Specifies πᵢ ∀ i ∈ {HPC, LPC}
    - Shaft power demand, Tt41 
    - DeNOx target
    - SCR paramters -> CPSI, w, l

"""
function NPSS_TShaft_input(NPSS, alt_in, MN_in, 
    SHP_dmd, Tt41, 
    LPC_PR, HPC_PR,
    cpsi, w, lcat, deNOx, first, LHV)

    input_string = "111 "*
    "first = $first;" *
    "Eng.Amb.alt_in = $(alt_in/ft_to_m);" *
    "Eng.Amb.MN_in  = $MN_in ;" *
   
    "Eng.ShP.HPX = $(SHP_dmd/hp_to_W) ;" * #convert to hp before sending to NPSS
    "Tt41 = $Tt41 ;" *
    "deNOx_target = $deNOx ;"*

    "Eng.FusEng.LHV = $(LHV*429.923); "*
    "Eng.CmpL.PRdes = $(LPC_PR) ;" *
    "Eng.CmpH.PRdes = $(HPC_PR) ;" *

    "Eng.PCEC.l = $lcat ;"*
    "Eng.PCEC.w = $w ;"*
    "Eng.PCEC.cpsi = $cpsi ;"*
    "\n"

    write(NPSS, input_string)

    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        ηtherm = parse(Float64, out[2] )
        mdotf  = parse(Float64, out[3] )
        BSFC   = parse(Float64, out[4] )
        deNOx  = parse(Float64, out[5] )
        mcat   = parse(Float64, out[6] )
        EINOx1 = parse(Float64, out[7] )
        EINOx2 = EINOx1*(1-deNOx)
        mdot   = parse(Float64, out[8] )
        Tt3    = parse(Float64, out[9] )
        OPR    = parse(Float64, out[10])
        Wc3    = parse(Float64, out[11])
        LHV    = parse(Float64, out[12])
        W_in    = parse(Float64, out[13])
    end
    # println(W_in, "lbm/s; If m/mdot^1.2 = 38.9 -> Wtshaft = ", (W_in/2.205)^1.2 * 38.9 * gee) # 38.9 kg/(kg/s) from Hall et al. https://arc.aiaa.org/doi/pdf/10.2514/6.2018-3973
    return NPSS_success, ηtherm, mdotf, BSFC, deNOx, mcat, EINOx1, EINOx2, mdot, Tt3, OPR, Wc3

end

"""
OffdesMode

Writes an input file for NPSS Turboshaft model in off-des conditions

    - Abmient altitude [in m]
    - Flight mach number
    - Tt41 
    - the flag first sets whether NPSS should restart or use previous value

"""
function NPSS_TShaft_run(NPSS, alt_in, MN_in, 
                            Tt41, N2_dmd, first)

    write(NPSS, "222 Tt41=$Tt41; N2_dmd = $N2_dmd; Eng.Amb.alt_in=$(alt_in/ft_to_m); Eng.Amb.MN_in = $MN_in;first=$first;\n")
    
    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        ShP    = parse(Float64, out[2] ) #NPSS returns in W
        ηtherm = parse(Float64, out[3] )
        mdotf  = parse(Float64, out[4] )
        BSFC   = parse(Float64, out[5] )
        deNOx  = parse(Float64, out[6] )
        EGT    = parse(Float64, out[7] )
        Tt3    = parse(Float64, out[8] )
        W3     = parse(Float64, out[9] )
        EINOx1 = parse(Float64, out[10] )
        EINOx2 = EINOx1*(1-deNOx)
        FAR    = parse(Float64, out[11])
    end
    
    return NPSS_success, ShP, ηtherm, mdotf, BSFC, deNOx, EGT, Tt3, W3, EINOx1, EINOx2, FAR
end

function NPSS_TShaft_run2(NPSS, alt_in, MN_in, 
    SHP_dmd, N2_dmd, first)

write(NPSS, "333 Eng.ShP.HPX=$(SHP_dmd/hp_to_W); N2_dmd = $N2_dmd; Eng.Amb.alt_in=$(alt_in/ft_to_m); Eng.Amb.MN_in = $MN_in; first=$first;\n")

out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
NPSS_success = parse(Float64, out[1] )

if length(out) > 1
ShP    = parse(Float64, out[2] )
ηtherm = parse(Float64, out[3] )
mdotf  = parse(Float64, out[4] )
BSFC   = parse(Float64, out[5] )
deNOx  = parse(Float64, out[6] )
end

return NPSS_success, ShP, ηtherm, mdotf, BSFC, deNOx
end
"""
Fan - Design mode:

Writes the desired conditions to the external NPSS process
"""
function runNPSS_Fan(NPSS::Base.Process, alt_in::Float64, MN_in::Float64, Fn::Float64,
                    Kinl::Float64, Φinl::Float64, 
                    πfan::Float64, first)
    input_string = "111 "*
                "first = $first;" *
                "DuctedFan.Amb.alt_in = $(alt_in/ft_to_m);" *
                "DuctedFan.Amb.MN_in  = $MN_in ;" *
                "Fn_target = $Fn ;" *
                "DuctedFan.InEng.Kinl = $Kinl;" *
                "DuctedFan.Phiinl     = $Φinl;" *
                "DuctedFan.Fan.PRdes = $πfan;" *
                "\n"

    write(NPSS, input_string)

    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        Dfan = parse(Float64, out[2] )
        Power  = parse(Float64, out[3] )
        Torque = parse(Float64, out[4] )
        Nmech  = parse(Float64, out[5] )
        Mtip   = parse(Float64, out[6] )
        ηprop  = parse(Float64, out[7] )
        ηDF    = parse(Float64, out[8] )
        ANoz   = parse(Float64, out[9] )
    end

    return NPSS_success, Dfan, Power, Torque, Nmech, Mtip, ηprop, ηDF, ANoz
end

"""
Fan - Off-design mode:

Writes the desired conditions to the external NPSS process
"""
function runNPSS_Fan(NPSS::Base.Process, alt_in::Float64, MN_in::Float64, Pin::Float64,
                    Kinl::Float64, Φinl::Float64, first)
    input_string = "222 "*
                "first = $first;" *
                "DuctedFan.Amb.alt_in = $(alt_in/ft_to_m);" *
                "DuctedFan.Amb.MN_in  = $MN_in ;" *
                "ShP_input = $(-Pin/745.7) ;" *
                "DuctedFan.InEng.Kinl = $Kinl;" *
                "DuctedFan.Phiinl     = $Φinl;" *
                "\n"

    write(NPSS, input_string)

    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        Fn_N = parse(Float64, out[2] )
        Power  = parse(Float64, out[3] )
        Torque = parse(Float64, out[4] )
        Nmech  = parse(Float64, out[5] )
        Mtip   = parse(Float64, out[6] )
        ηprop  = parse(Float64, out[7] )
        ηDF    = parse(Float64, out[8] )
    end

    return NPSS_success, Fn_N, Power, Torque, Nmech, Mtip, ηprop, ηDF
end
function runNPSS_Fan2(NPSS::Base.Process, alt_in::Float64, MN_in::Float64, Fn::Float64,
    Kinl::Float64, Φinl::Float64, first)

    input_string = "333 "*
"first = $first;" *
"DuctedFan.Amb.alt_in = $(alt_in/ft_to_m);" *
"DuctedFan.Amb.MN_in  = $MN_in ;" *
"Fn_target = $Fn ;" *
"DuctedFan.InEng.Kinl = $Kinl;" *
"DuctedFan.Phiinl     = $Φinl;" *
"\n"

write(NPSS, input_string)

out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
NPSS_success = parse(Float64, out[1] )

if length(out) > 1
Fn_N = parse(Float64, out[2] )
Power  = parse(Float64, out[3] )
Torque = parse(Float64, out[4] )
Nmech  = parse(Float64, out[5] )
Mtip   = parse(Float64, out[6] )
ηprop  = parse(Float64, out[7] )
ηDF    = parse(Float64, out[8] )
end

return NPSS_success, Fn_N, Power, Torque, Nmech, Mtip, ηprop, ηDF
end