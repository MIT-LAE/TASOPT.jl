"""
`weight_buildup` prints out the weight build up for the aircraft
"""
function weight_buildup(ac::aircraft; io=stdout)
    parg = ac.parg
    pari = ac.pari
    fuselage = ac.fuselage
    Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
    Whpesys = parg[igWMTO] * fuselage.HPE_sys.W
    Wlgnose = parg[igWMTO] * parg[igflgnose]
    Wlgmain = parg[igWMTO] * parg[igflgmain]
    Wtotadd = Whpesys + Wlgnose + Wlgmain

    Wbox    = parg[igWweb] + parg[igWcap]
    Wflap   = Wbox * parg[igfflap]
    Wslat   = Wbox * parg[igfslat]
    Waile   = Wbox * parg[igfaile]
    Wlete   = Wbox * parg[igflete]
    Wribs   = Wbox * parg[igfribs]
    Wspoi   = Wbox * parg[igfspoi]
    Wwatt   = Wbox * parg[igfwatt]
    Wwing   = Wbox + Wflap + Wslat + Waile + Wlete + Wribs + Wspoi + Wwatt


    printstyled(io, "Weight build-up:\n -------------- \n", color=:bold )
    @printf(io, "Wempty  + %10.1f N (%8.1f lb)\n", Wempty, Wempty/lbf_to_N)
    @printf(io, "Wpay    + %10.1f N (%8.1f lb)\n", parg[igWpay], parg[igWpay]/lbf_to_N)
    @printf(io, "Wfuel   + %10.1f N (%8.1f lb)\n", parg[igWfuel], parg[igWfuel]/lbf_to_N)
    @printf(io, "--------------------\n")
    printstyled(io, @sprintf("WMTO    = %10.1f N (%8.1f lb)\n\n",
                         parg[igWMTO], parg[igWMTO]/lbf_to_N); color=:bold)

    @printf(io,"Wfuse   + %10.1f N (%8.1f lb)\n", fuselage.weight, fuselage.weight/lbf_to_N)
    @printf(io,"Wwing   + %10.1f N (%8.1f lb)\n", parg[igWwing ], parg[igWwing ]/lbf_to_N)
    @printf(io,"Wvtail  + %10.1f N (%8.1f lb)\n", parg[igWvtail], parg[igWvtail]/lbf_to_N)
    @printf(io,"Whtail  + %10.1f N (%8.1f lb)\n", parg[igWhtail], parg[igWhtail]/lbf_to_N)
    @printf(io,"Wtesys  + %10.1f N (%8.1f lb)\n", parg[igWtesys], parg[igWtesys]/lbf_to_N)
    @printf(io,"Wftank  + %10.1f N (%8.1f lb)\n", parg[igWftank], parg[igWftank]/lbf_to_N)
    @printf(io,"Wadd    + %10.1f N (%8.1f lb)\n", Wtotadd, Wtotadd/lbf_to_N)
    @printf(io,"--------------------\n")
    printstyled(io, @sprintf("Wempty  = %10.1f N (%8.1f lb)\n\n", 
    fuselage.weight + parg[igWwing]+ parg[igWvtail] + parg[igWhtail] + 
    parg[igWtesys] + +parg[igWftank] + Wtotadd, 
    (fuselage.weight + parg[igWwing]+ parg[igWvtail] + parg[igWhtail] + 
    parg[igWtesys] + +parg[igWftank] + Wtotadd)/lbf_to_N); color=:bold)

    @printf(io,"Wcap    + %10.1f N (%8.1f lb)\n", parg[igWcap], parg[igWcap]/lbf_to_N)
    @printf(io,"Wweb    + %10.1f N (%8.1f lb)\n", parg[igWweb], parg[igWweb]/lbf_to_N)
    @printf(io,"Wflap   + %10.1f N (%8.1f lb)\n", Wflap, Wflap/lbf_to_N)
    @printf(io,"Wslat   + %10.1f N (%8.1f lb)\n", Wslat, Wslat/lbf_to_N)
    @printf(io,"Waile   + %10.1f N (%8.1f lb)\n", Waile, Waile/lbf_to_N)
    @printf(io,"Wlete   + %10.1f N (%8.1f lb)\n", Wlete, Wlete/lbf_to_N)
    @printf(io,"Wribs   + %10.1f N (%8.1f lb)\n", Wribs, Wribs/lbf_to_N)
    @printf(io,"Wspoi   + %10.1f N (%8.1f lb)\n", Wspoi, Wspoi/lbf_to_N)
    @printf(io,"Wwatt   + %10.1f N (%8.1f lb)\n", Wwatt, Wwatt/lbf_to_N)
    @printf(io,"--------------------\n")
    printstyled(io, @sprintf("Wwing  = %10.1f N (%8.1f lb)\n\n", Wwing, Wwing/lbf_to_N); color=:bold)

    @printf(io,"Wtshaft + %10.1f N (%8.1f lb) × %d\n", parg[igWtshaft], parg[igWtshaft]/lbf_to_N, parpt[ipt_nTshaft])
    @printf(io,"Wcat    + %10.1f N (%8.1f lb) × %d\n", parg[igWcat   ], parg[igWcat   ]/lbf_to_N, parpt[ipt_nTshaft])
    @printf(io,"Waftfan + %10.1f N (%8.1f lb) × %d\n", parg[igWaftfan], parg[igWaftfan]/lbf_to_N, 2.0) 
    @printf(io,"WaftGB  + %10.1f N (%8.1f lb) × %d\n", parg[igWaftfanGB ], parg[igWaftfanGB ]/lbf_to_N, parpt[ipt_nTshaft]) 
    @printf(io,"Wgen    + %10.1f N (%8.1f lb) × %d\n", parg[igWgen   ], parg[igWgen   ]/lbf_to_N, parpt[ipt_ngen]) 
    @printf(io,"Wrect   + %10.1f N (%8.1f lb) × %d\n", parg[igWrect  ], parg[igWrect  ]/lbf_to_N, parpt[ipt_ngen]) 
    @printf(io,"Wcables + %10.1f N (%8.1f lb) ----\n", parg[igWcables], parg[igWcables]/lbf_to_N) 
    @printf(io,"Winv    + %10.1f N (%8.1f lb) × %d\n", parg[igWinv   ], parg[igWinv   ]/lbf_to_N, parpt[ipt_nfan]) 
    @printf(io,"Wmot    + %10.1f N (%8.1f lb) × %d\n", parg[igWmot   ], parg[igWmot   ]/lbf_to_N, parpt[ipt_nfan]) 
    @printf(io,"Wfan    + %10.1f N (%8.1f lb) × %d\n", parg[igWfan   ], parg[igWfan   ]/lbf_to_N, parpt[ipt_nfan]) 
    @printf(io,"WfanGB  + %10.1f N (%8.1f lb) × %d\n", parg[igWfanGB ], parg[igWfanGB ]/lbf_to_N, parpt[ipt_nfan]) 
    @printf(io,"Wtms    + %10.1f N (%8.1f lb) ----\n", parg[igWtms   ], parg[igWtms   ]/lbf_to_N) 
    @printf(io,"--------------------\n")
    printstyled(io,@sprintf("Wtesys  = %10.1f N (%8.1f lb)\n\n",
     parg[igWtesys], parg[igWtesys]/lbf_to_N ), color = :bold)

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
    @views para = ac.para[:,:,1]
    printstyled(io, "Aerodynamics:\n -------------- \n", color=:bold)

    @printf(io, "Ref.Area= %6.5f m²\n", parg[igS])
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
    @printf(io, "CL     × Sref = %6.5f m²\n", parg[igS]*para[iaCL, ipcruise1])
    @printf(io, "CD     × Sref = %6.5f m²\n", parg[igS]*para[iaCD, ipcruise1])
    @printf(io, "CDfuse × Sref = %6.5f m²\n", parg[igS]*para[iaCDfuse, ipcruise1])
    @printf(io, "CDi    × Sref = %6.5f m²\n", parg[igS]*para[iaCDi, ipcruise1])
    @printf(io, "CDwing × Sref = %6.5f m²\n", parg[igS]*para[iaCDwing, ipcruise1])
    @printf(io, "CDhtail× Sref = %6.5f m²\n", parg[igS]*para[iaCDhtail, ipcruise1])
    @printf(io, "CDvtail× Sref = %6.5f m²\n", parg[igS]*para[iaCDvtail, ipcruise1])
    @printf(io, "CDnace × Sref = %6.5f m²\n", parg[igS]*para[iaCDnace, ipcruise1])
    @printf(io, "CDBLIf × Sref = %6.5f m²\n", parg[igS]*para[iadCDBLIf, ipcruise1])
    @printf(io, "CDBLIw × Sref = %6.5f m²\n\n", parg[igS]*para[iadCDBLIw, ipcruise1])

end

"""
`geometry` prints out a numerical description of the aircraft layout
"""
function geometry(ac::aircraft; io = stdout)
    parg = ac.parg
    fuselage = ac.fuselage
    printstyled(io, "Fuselage Layout:\n -------------- \n", color=:bold )
    @printf(io, "xnose     = %5.1f m (%8.1f ft)\n", fuselage.layout.x_nose , fuselage.layout.x_nose/ft_to_m)
    @printf(io, "xend      = %5.1f m (%8.1f ft)\n", fuselage.layout.x_end , fuselage.layout.x_end/ft_to_m)
    @printf(io, "xwing     = %5.1f m (%8.1f ft)\n", parg[igxwing  ] , parg[igxwing   ]/ft_to_m)
    @printf(io, "xhtail    = %5.1f m (%8.1f ft)\n", parg[igxhtail ] , parg[igxhtail  ]/ft_to_m)
    @printf(io, "xvtail    = %5.1f m (%8.1f ft)\n", parg[igxvtail ] , parg[igxvtail  ]/ft_to_m)
    @printf(io, "xblend1   = %5.1f m (%8.1f ft)\n", fuselage.layout.x_start_cylinder , fuselage.layout.x_start_cylinder/ft_to_m)
    @printf(io, "xblend2   = %5.1f m (%8.1f ft)\n", fuselage.layout.x_end_cylinder , fuselage.layout.x_end_cylinder/ft_to_m)
    @printf(io, "xshell1   = %5.1f m (%8.1f ft)\n", fuselage.layout.x_pressure_shell_fwd , fuselage.layout.x_pressure_shell_fwd/ft_to_m)
    @printf(io, "xshell2   = %5.1f m (%8.1f ft)\n", fuselage.layout.x_pressure_shell_aft , fuselage.layout.x_pressure_shell_aft/ft_to_m)
    @printf(io, "xhbend    = %5.1f m (%8.1f ft)\n", fuselage.bendingmaterial_h.weight.x , fuselage.bendingmaterial_h.weight.x/ft_to_m)
    @printf(io, "xvbend    = %5.1f m (%8.1f ft)\n", fuselage.bendingmaterial_v.weight.x , fuselage.bendingmaterial_v.weight.x/ft_to_m)
    @printf(io, "xwbox     = %5.1f m (%8.1f ft)\n", parg[igxwbox  ] , parg[igxwbox   ]/ft_to_m)
    @printf(io, "xhbox     = %5.1f m (%8.1f ft)\n", parg[igxhbox  ] , parg[igxhbox   ]/ft_to_m)
    @printf(io, "xvbox     = %5.1f m (%8.1f ft)\n", parg[igxvbox  ] , parg[igxvbox   ]/ft_to_m)
    @printf(io, "xtshaft   = %5.1f m (%8.1f ft)\n", parg[igxtshaft] , parg[igxtshaft ]/ft_to_m)
    @printf(io, "xgen      = %5.1f m (%8.1f ft)\n", parg[igxgen   ] , parg[igxgen    ]/ft_to_m)
    @printf(io, "xcat      = %5.1f m (%8.1f ft)\n", parg[igxcat   ] , parg[igxcat    ]/ft_to_m)
    @printf(io, "xftank    = %5.1f m (%8.1f ft)\n", parg[igxftank ] , parg[igxftank  ]/ft_to_m)
    @printf(io, "xftankaft = %5.1f m (%8.1f ft)\n", parg[igxftankaft ] , parg[igxftankaft  ]/ft_to_m)
    
    @printf(io, "\nRfuse  = %5.1f m (%8.1f ft)\n", fuselage.layout.cross_section.radius , fuselage.layout.cross_section.radius/ft_to_m)

    
    SMfwd = (parg[igxNP] - parg[igxCGfwd])/parg[igcma]
    SMaft = (parg[igxNP] - parg[igxCGaft])/parg[igcma]
    printstyled(io, "\nStability:\n -------------- \n", color=:bold )
    @printf(io, "xNP     = %5.1f m (%8.1f ft)\n", parg[igxNP ] , parg[igxNP ]/ft_to_m)
    @printf(io, "xCGfwd  = %5.1f m (%8.1f ft)\n", parg[igxCGfwd ] , parg[igxCGfwd ]/ft_to_m)
    @printf(io, "xCGaft  = %5.1f m (%8.1f ft)\n", parg[igxCGaft ] , parg[igxCGaft ]/ft_to_m)
    @printf(io, "xSMfwd  = %5.4f\n", SMfwd)
    @printf(io, "xSMaft  = %5.4f\n", SMaft)

    
    printstyled(io, "\nWing Layout:\n -------------- \n", color=:bold )
    @printf(io, "AR      = %5.3f \n" , parg[igAR     ])
    @printf(io, "sweep   = %5.3f \n" , parg[igsweep  ])
    @printf(io, "lambdas = %5.3f \n" , parg[iglambdas])
    @printf(io, "lambdat = %5.3f \n" , parg[iglambdat]) 
    co = parg[igco]
    cs = parg[igco]*parg[iglambdas]
    ct = parg[igco]*parg[iglambdat]

    @printf(io, "co      = %5.1f m (%8.1f ft)\n" , co, co / ft_to_m )
    @printf(io, "cs      = %5.1f m (%8.1f ft)\n" , cs, cs / ft_to_m )
    @printf(io, "ct      = %5.1f m (%8.1f ft)\n" , ct, ct / ft_to_m )
    @printf(io, "bo      = %5.1f m (%8.1f ft)\n" , parg[igbo], parg[igbo]/ft_to_m   )
    @printf(io, "bs      = %5.1f m (%8.1f ft)\n" , parg[igbs], parg[igbs]/ft_to_m   )
    @printf(io, "b       = %5.1f m (%8.1f ft)\n" , parg[igb ], parg[igb ]/ft_to_m   )
    @printf(io, "S       = %5.1f m²(%8.1f ft²)\n" , parg[igS ], parg[igS ]/ft_to_m^2 )


end