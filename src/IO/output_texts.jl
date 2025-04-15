"""
`weight_buildup` prints out the weight build up for the aircraft
"""
function weight_buildup(ac::aircraft; io=stdout)
    parg = ac.parg
    fuselage = ac.fuselage
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    landing_gear = ac.landing_gear
    Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
    Whpesys = parg[igWMTO] * fuselage.HPE_sys.W
    Wlgnose = landing_gear.nose_gear.weight.W
    Wlgmain = landing_gear.main_gear.weight.W
    Wtotadd = Whpesys + Wlgnose + Wlgmain

    Wbox    = wing.inboard.webs.weight.W + wing.inboard.caps.weight.W
    Wflap   = Wbox * wing.weight_frac_flap
    Wslat   = Wbox * wing.weight_frac_slat
    Waile   = Wbox * wing.weight_frac_ailerons
    Wlete   = Wbox * wing.weight_frac_leading_trailing_edge
    Wribs   = Wbox * wing.weight_frac_ribs
    Wspoi   = Wbox * wing.weight_frac_spoilers
    Wwatt   = Wbox * wing.weight_frac_attachments
    Wwing   = Wbox + Wflap + Wslat + Waile + Wlete + Wribs + Wspoi + Wwatt


    printstyled(io, "Weight build-up:\n -------------- \n", color=:bold )
    @printf(io, "Wempty  + %10.1f N (%8.1f lb)\n", Wempty, Wempty/lbf_to_N)
    @printf(io, "Wpay    + %10.1f N (%8.1f lb)\n", parg[igWpay], parg[igWpay]/lbf_to_N)
    @printf(io, "Wfuel   + %10.1f N (%8.1f lb)\n", parg[igWfuel], parg[igWfuel]/lbf_to_N)
    @printf(io, "--------------------\n")
    printstyled(io, @sprintf("WMTO    = %10.1f N (%8.1f lb)\n\n",
                         parg[igWMTO], parg[igWMTO]/lbf_to_N); color=:bold)

    @printf(io,"Wfuse   + %10.1f N (%8.1f lb)\n", fuselage.weight, fuselage.weight/lbf_to_N)
    @printf(io,"Wwing   + %10.1f N (%8.1f lb)\n", wing.weight, wing.weight/lbf_to_N)
    @printf(io,"Wvtail  + %10.1f N (%8.1f lb)\n", vtail.weight, vtail.weight/lbf_to_N)
    @printf(io,"Whtail  + %10.1f N (%8.1f lb)\n", htail.weight, htail.weight/lbf_to_N)
    @printf(io,"Wtesys  + %10.1f N (%8.1f lb)\n", parg[igWtesys], parg[igWtesys]/lbf_to_N)
    @printf(io,"Wftank  + %10.1f N (%8.1f lb)\n", parg[igWftank], parg[igWftank]/lbf_to_N)
    @printf(io,"Wadd    + %10.1f N (%8.1f lb)\n", Wtotadd, Wtotadd/lbf_to_N)
    @printf(io,"--------------------\n")
    printstyled(io, @sprintf("Wempty  = %10.1f N (%8.1f lb)\n\n", 
    fuselage.weight + wing.weight + vtail.weight + htail.weight + 
    parg[igWtesys] + +parg[igWftank] + Wtotadd, 
    (fuselage.weight + wing.weight + vtail.weight + htail.weight + 
    parg[igWtesys] + +parg[igWftank] + Wtotadd)/lbf_to_N); color=:bold)

    @printf(io,"Wcap    + %10.1f N (%8.1f lb)\n", wing.inboard.caps.weight.W, wing.inboard.caps.weight.W/lbf_to_N)
    @printf(io,"Wweb    + %10.1f N (%8.1f lb)\n", wing.inboard.webs.weight.W , wing.inboard.webs.weight.W /lbf_to_N)
    @printf(io,"Wflap   + %10.1f N (%8.1f lb)\n", Wflap, Wflap/lbf_to_N)
    @printf(io,"Wslat   + %10.1f N (%8.1f lb)\n", Wslat, Wslat/lbf_to_N)
    @printf(io,"Waile   + %10.1f N (%8.1f lb)\n", Waile, Waile/lbf_to_N)
    @printf(io,"Wlete   + %10.1f N (%8.1f lb)\n", Wlete, Wlete/lbf_to_N)
    @printf(io,"Wribs   + %10.1f N (%8.1f lb)\n", Wribs, Wribs/lbf_to_N)
    @printf(io,"Wspoi   + %10.1f N (%8.1f lb)\n", Wspoi, Wspoi/lbf_to_N)
    @printf(io,"Wwatt   + %10.1f N (%8.1f lb)\n", Wwatt, Wwatt/lbf_to_N)
    @printf(io,"--------------------\n")
    printstyled(io, @sprintf("Wwing  = %10.1f N (%8.1f lb)\n\n", Wwing, Wwing/lbf_to_N); color=:bold)

    @printf(io,"Wftnkins  = %10.1f N (%8.1f lb)\n", parg[igWinsftank ], parg[igWinsftank ]/lbf_to_N) 
    @printf(io,"Wftank    = %10.1f N (%8.1f lb)\n\n", parg[igWftank    ], parg[igWftank    ]/lbf_to_N) 
    @printf(io,"lftank    = %10.1f m (%8.1f ft)\n"  , parg[iglftank    ], parg[iglftank    ]/ft_to_m) 
    @printf(io,"Rftank    = %10.1f m (%8.1f ft)\n"  , parg[igRftank    ], parg[igRftank    ]/ft_to_m) 

    @printf(io,"ηtank     = %3.1f %% \n\n", parg[igWfuel]/(parg[igWfuel] + parg[igWftank])*100)
end
"""
    aero(parg, para; io = stdout)

`aero` returns a summary of all aerodynamic properties 
of the aircraft
"""
function aero(ac::aircraft; io = stdout)
    parg = ac.parg
    wing = ac.wing
    @views para = ac.para[:,:,1]
    printstyled(io, "Aerodynamics:\n -------------- \n", color=:bold)

    @printf(io, "Ref.Area= %6.5f m²\n", wing.layout.S)
    @printf(io, "L/D     = %6.5f\n", para[iaCL, ipcruise1]/ para[iaCD, ipcruise1])
    @printf(io, "CL      = %6.5f\n", para[iaCL, ipcruise1])
    @printf(io, "CD      = %6.5f\n", para[iaCD, ipcruise1])
    @printf(io, "CDfuse  = %6.5f\n", para[iaCDfuse, ipcruise1])
    @printf(io, "CDi     = %6.5f\n", para[iaCDi, ipcruise1])
    @printf(io, "CDwing  = %6.5f\n", para[iaCDwing, ipcruise1])
    @printf(io, "CDhtail = %6.5f\n", para[iaCDhtail, ipcruise1])
    @printf(io, "CDvtail = %6.5f\n", para[iaCDvtail, ipcruise1])
    @printf(io, "CDnace  = %6.5f\n", para[iaCDnace, ipcruise1])
    @printf(io, "CDBLIf  = %6.5f\n", para[iadCDBLIf, ipcruise1])
    @printf(io, "CDBLIw  = %6.5f\n", para[iadCDBLIw, ipcruise1])

    printstyled(io, "\nDrag Areas = CD × Sref:\n", color=:bold)
    @printf(io, "CL     × Sref = %6.5f m²\n", wing.layout.S*para[iaCL, ipcruise1])
    @printf(io, "CD     × Sref = %6.5f m²\n", wing.layout.S*para[iaCD, ipcruise1])
    @printf(io, "CDfuse × Sref = %6.5f m²\n", wing.layout.S*para[iaCDfuse, ipcruise1])
    @printf(io, "CDi    × Sref = %6.5f m²\n", wing.layout.S*para[iaCDi, ipcruise1])
    @printf(io, "CDwing × Sref = %6.5f m²\n", wing.layout.S*para[iaCDwing, ipcruise1])
    @printf(io, "CDhtail× Sref = %6.5f m²\n", wing.layout.S*para[iaCDhtail, ipcruise1])
    @printf(io, "CDvtail× Sref = %6.5f m²\n", wing.layout.S*para[iaCDvtail, ipcruise1])
    @printf(io, "CDnace × Sref = %6.5f m²\n", wing.layout.S*para[iaCDnace, ipcruise1])
    @printf(io, "CDBLIf × Sref = %6.5f m²\n", wing.layout.S*para[iadCDBLIf, ipcruise1])
    @printf(io, "CDBLIw × Sref = %6.5f m²\n\n", wing.layout.S*para[iadCDBLIw, ipcruise1])

end

"""
`geometry` prints out a numerical description of the aircraft layout
"""
function geometry(ac::aircraft; io = stdout)
    parg = ac.parg
    fuselage = ac.fuselage
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    printstyled(io, "Fuselage Layout:\n -------------- \n", color=:bold )
    @printf(io, "xnose     = %5.1f m (%8.1f ft)\n", fuselage.layout.x_nose , fuselage.layout.x_nose/ft_to_m)
    @printf(io, "xend      = %5.1f m (%8.1f ft)\n", fuselage.layout.x_end , fuselage.layout.x_end/ft_to_m)
    @printf(io, "xwing     = %5.1f m (%8.1f ft)\n", wing.layout.x , wing.layout.x/ft_to_m)
    @printf(io, "xhtail    = %5.1f m (%8.1f ft)\n", htail.layout.x , htail.layout.x/ft_to_m)
    @printf(io, "xvtail    = %5.1f m (%8.1f ft)\n", vtail.layout.x , vtail.layout.x/ft_to_m)
    @printf(io, "xblend1   = %5.1f m (%8.1f ft)\n", fuselage.layout.x_start_cylinder , fuselage.layout.x_start_cylinder/ft_to_m)
    @printf(io, "xblend2   = %5.1f m (%8.1f ft)\n", fuselage.layout.x_end_cylinder , fuselage.layout.x_end_cylinder/ft_to_m)
    @printf(io, "xshell1   = %5.1f m (%8.1f ft)\n", fuselage.layout.x_pressure_shell_fwd , fuselage.layout.x_pressure_shell_fwd/ft_to_m)
    @printf(io, "xshell2   = %5.1f m (%8.1f ft)\n", fuselage.layout.x_pressure_shell_aft , fuselage.layout.x_pressure_shell_aft/ft_to_m)
    @printf(io, "xhbend    = %5.1f m (%8.1f ft)\n", fuselage.bendingmaterial_h.weight.x , fuselage.bendingmaterial_h.weight.x/ft_to_m)
    @printf(io, "xvbend    = %5.1f m (%8.1f ft)\n", fuselage.bendingmaterial_v.weight.x , fuselage.bendingmaterial_v.weight.x/ft_to_m)
    @printf(io, "xwbox     = %5.1f m (%8.1f ft)\n", wing.layout.box_x , wing.layout.box_x/ft_to_m)
    @printf(io, "xhbox     = %5.1f m (%8.1f ft)\n", htail.layout.box_x , htail.layout.box_x/ft_to_m)
    @printf(io, "xvbox     = %5.1f m (%8.1f ft)\n", vtail.layout.box_x , vtail.layout.box_x/ft_to_m)
    @printf(io, "xtshaft   = %5.1f m (%8.1f ft)\n", parg[igxtshaft] , parg[igxtshaft ]/ft_to_m)
    @printf(io, "xgen      = %5.1f m (%8.1f ft)\n", parg[igxgen   ] , parg[igxgen    ]/ft_to_m)
    @printf(io, "xcat      = %5.1f m (%8.1f ft)\n", parg[igxcat   ] , parg[igxcat    ]/ft_to_m)
    @printf(io, "xftank    = %5.1f m (%8.1f ft)\n", parg[igxftank ] , parg[igxftank  ]/ft_to_m)
    @printf(io, "xftankaft = %5.1f m (%8.1f ft)\n", parg[igxftankaft ] , parg[igxftankaft  ]/ft_to_m)
    
    @printf(io, "\nRfuse  = %5.1f m (%8.1f ft)\n", fuselage.layout.cross_section.radius , fuselage.layout.cross_section.radius/ft_to_m)

    
    SMfwd = (parg[igxNP] - parg[igxCGfwd])/wing.mean_aero_chord
    SMaft = (parg[igxNP] - parg[igxCGaft])/wing.mean_aero_chord
    printstyled(io, "\nStability:\n -------------- \n", color=:bold )
    @printf(io, "xNP     = %5.1f m (%8.1f ft)\n", parg[igxNP ] , parg[igxNP ]/ft_to_m)
    @printf(io, "xCGfwd  = %5.1f m (%8.1f ft)\n", parg[igxCGfwd ] , parg[igxCGfwd ]/ft_to_m)
    @printf(io, "xCGaft  = %5.1f m (%8.1f ft)\n", parg[igxCGaft ] , parg[igxCGaft ]/ft_to_m)
    @printf(io, "xSMfwd  = %5.4f\n", SMfwd)
    @printf(io, "xSMaft  = %5.4f\n", SMaft)

    printstyled(io, "\nWing Layout:\n -------------- \n", color=:bold )
    @printf(io, "AR      = %5.3f \n" , wing.layout.AR)
    @printf(io, "sweep   = %5.3f \n" , wing.layout.sweep)
    @printf(io, "lambdas = %5.3f \n" , wing.inboard.λ)
    @printf(io, "lambdat = %5.3f \n" , wing.outboard.λ) 
    co = wing.layout.root_chord
    cs = wing.layout.root_chord*wing.inboard.λ 
    ct = wing.layout.root_chord*wing.outboard.λ 

    @printf(io, "co      = %5.1f m (%8.1f ft)\n" , co, co / ft_to_m )
    @printf(io, "cs      = %5.1f m (%8.1f ft)\n" , cs, cs / ft_to_m )
    @printf(io, "ct      = %5.1f m (%8.1f ft)\n" , ct, ct / ft_to_m )
    @printf(io, "bo      = %5.1f m (%8.1f ft)\n" , wing.layout.root_span, wing.layout.root_span/ft_to_m   )
    @printf(io, "bs      = %5.1f m (%8.1f ft)\n" , wing.layout.break_span, wing.layout.break_span/ft_to_m   )
    @printf(io, "b       = %5.1f m (%8.1f ft)\n" , wing.layout.span, wing.layout.span/ft_to_m   )
    @printf(io, "S       = %5.1f m²(%8.1f ft²)\n" , wing.layout.S, wing.layout.S/ft_to_m^2 )


end