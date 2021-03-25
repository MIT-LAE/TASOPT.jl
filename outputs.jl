"""
`weight_buildup` prints out the weight build up for the aircraft
"""
function weight_buildup(parg; io=stdout)

    Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
    Whpesys = parg[igWMTO] * parg[igfhpesys]
    Wlgnose = parg[igWMTO] * parg[igflgnose]
    Wlgmain = parg[igWMTO] * parg[igflgmain]
    Wtotadd = Whpesys + Wlgnose + Wlgmain

    printstyled(io, "Weight build-up:\n -------------- \n", color=:bold )
    @printf(io, "Wempty  + %10.1f N (%8.1f lb)\n", Wempty, Wempty/lbf_to_N)
    @printf(io, "Wpay    + %10.1f N (%8.1f lb)\n", parg[igWpay], parg[igWpay]/lbf_to_N)
    @printf(io, "Wfuel   + %10.1f N (%8.1f lb)\n", parg[igWfuel], parg[igWfuel]/lbf_to_N)
    @printf(io, "--------------------\n")
    printstyled(io, @sprintf("WMTO    = %10.1f N (%8.1f lb)\n\n",
                         parg[igWMTO], parg[igWMTO]/lbf_to_N); color=:bold)

    @printf(io,"Wfuse   + %10.1f N (%8.1f lb)\n", parg[igWfuse ], parg[igWfuse ]/lbf_to_N)
    @printf(io,"Wwing   + %10.1f N (%8.1f lb)\n", parg[igWwing ], parg[igWwing ]/lbf_to_N)
    @printf(io,"Wvtail  + %10.1f N (%8.1f lb)\n", parg[igWvtail], parg[igWvtail]/lbf_to_N)
    @printf(io,"Whtail  + %10.1f N (%8.1f lb)\n", parg[igWhtail], parg[igWhtail]/lbf_to_N)
    @printf(io,"Wtesys  + %10.1f N (%8.1f lb)\n", parg[igWtesys], parg[igWtesys]/lbf_to_N)
    @printf(io,"Wftank  + %10.1f N (%8.1f lb)\n", parg[igWftank], parg[igWftank]/lbf_to_N)
    @printf(io,"Wadd    + %10.1f N (%8.1f lb)\n", Wtotadd, Wtotadd/lbf_to_N)
    @printf(io,"--------------------\n")
    printstyled(io, @sprintf("Wempty  = %10.1f N (%8.1f lb)\n\n", 
    parg[igWfuse] + parg[igWwing]+ parg[igWvtail] + parg[igWhtail] + 
    parg[igWtesys] + +parg[igWftank] + Wtotadd, 
    (parg[igWfuse] + parg[igWwing]+ parg[igWvtail] + parg[igWhtail] + 
    parg[igWtesys] + +parg[igWftank] + Wtotadd)/lbf_to_N); color=:bold)

    @printf(io,"Wtshaft + %10.1f N × %d\n", parg[igWtshaft], parpt[ipt_nTshaft])
    @printf(io,"Wcat    + %10.1f N × %d\n", parg[igWcat   ], parpt[ipt_nTshaft])
    @printf(io,"Wgen    + %10.1f N × %d\n", parg[igWgen   ], parpt[ipt_ngen]) 
    @printf(io,"Winv    + %10.1f N × %d\n", parg[igWinv   ], parpt[ipt_nfan]) 
    @printf(io,"Wmot    + %10.1f N × %d\n", parg[igWmot   ], parpt[ipt_nfan]) 
    @printf(io,"Wfan    + %10.1f N × %d\n", parg[igWfan   ], parpt[ipt_nfan]) 
    @printf(io,"--------------------\n")
    printstyled(io,@sprintf("Wtesys  = %10.1f N (%8.1f lb)\n\n",
     parg[igWtesys], parg[igWtesys]/lbf_to_N ), color = :bold)

    @printf(io,"Wftnkins  = %10.1f N (%8.1f lb)\n", parg[igWinsftank ], parg[igWinsftank ]/lbf_to_N) 
    @printf(io,"Wftank    = %10.1f N (%8.1f lb)\n\n", parg[igWftank    ], parg[igWftank    ]/lbf_to_N) 
    @printf(io,"lftank    = %10.1f m (%8.1f ft)\n"  , parg[iglftank    ], parg[iglftank    ]/ft_to_m) 
    @printf(io,"Rftank    = %10.1f m (%8.1f ft)\n"  , parg[igRftank    ], parg[igRftank    ]/ft_to_m) 

    @printf(io,"ηtank   = %3.1f %% \n\n", parg[igWfuel]/(parg[igWfuel] + parg[igWftank])*100)
end

"""
`geometry` prints out the layout of the aircraft
"""
function geometry(parg; io = stdout)

    printstyled(io, "Fuselage Layout:\n -------------- \n", color=:bold )
    @printf(io, "xnose   = %5.1f m (%8.1f ft)\n", parg[igxnose  ] , parg[igxnose   ]/ft_to_m)
    @printf(io, "xend    = %5.1f m (%8.1f ft)\n", parg[igxend   ] , parg[igxend    ]/ft_to_m)
    @printf(io, "xwing   = %5.1f m (%8.1f ft)\n", parg[igxwing  ] , parg[igxwing   ]/ft_to_m)
    @printf(io, "xhtail  = %5.1f m (%8.1f ft)\n", parg[igxhtail ] , parg[igxhtail  ]/ft_to_m)
    @printf(io, "xvtail  = %5.1f m (%8.1f ft)\n", parg[igxvtail ] , parg[igxvtail  ]/ft_to_m)
    @printf(io, "xblend1 = %5.1f m (%8.1f ft)\n", parg[igxblend1] , parg[igxblend1 ]/ft_to_m)
    @printf(io, "xblend2 = %5.1f m (%8.1f ft)\n", parg[igxblend2] , parg[igxblend2 ]/ft_to_m)
    @printf(io, "xshell1 = %5.1f m (%8.1f ft)\n", parg[igxshell1] , parg[igxshell1 ]/ft_to_m)
    @printf(io, "xshell2 = %5.1f m (%8.1f ft)\n", parg[igxshell2] , parg[igxshell2 ]/ft_to_m)
    @printf(io, "xhbend  = %5.1f m (%8.1f ft)\n", parg[igxhbend ] , parg[igxhbend  ]/ft_to_m)
    @printf(io, "xvbend  = %5.1f m (%8.1f ft)\n", parg[igxvbend ] , parg[igxvbend  ]/ft_to_m)
    @printf(io, "xwbox   = %5.1f m (%8.1f ft)\n", parg[igxwbox  ] , parg[igxwbox   ]/ft_to_m)
    @printf(io, "xhbox   = %5.1f m (%8.1f ft)\n", parg[igxhbox  ] , parg[igxhbox   ]/ft_to_m)
    @printf(io, "xvbox   = %5.1f m (%8.1f ft)\n", parg[igxvbox  ] , parg[igxvbox   ]/ft_to_m)
    @printf(io, "xtshaft = %5.1f m (%8.1f ft)\n", parg[igxtshaft] , parg[igxtshaft ]/ft_to_m)
    @printf(io, "xgen    = %5.1f m (%8.1f ft)\n", parg[igxgen   ] , parg[igxgen    ]/ft_to_m)
    @printf(io, "xcat    = %5.1f m (%8.1f ft)\n", parg[igxcat   ] , parg[igxcat    ]/ft_to_m)
    @printf(io, "xftank  = %5.1f m (%8.1f ft)\n", parg[igxftank ] , parg[igxftank  ]/ft_to_m)

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