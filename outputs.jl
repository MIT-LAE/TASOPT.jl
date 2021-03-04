"""
`weight_buildup` prints out the weight build up for the aircraft
"""
function weight_buildup(parg)

    Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
    Whpesys = parg[igWMTO] * parg[igfhpesys]
    Wlgnose = parg[igWMTO] * parg[igflgnose]
    Wlgmain = parg[igWMTO] * parg[igflgmain]
    Wtotadd = Whpesys + Wlgnose + Wlgmain

    printstyled("Weight build-up:\n -------------- \n", color=:bold )
    @printf("Wempty  + %10.1f N (%8.1f lb)\n", Wempty, Wempty/lbf_to_N)
    @printf("Wpay    + %10.1f N (%8.1f lb)\n", parg[igWpay], parg[igWpay]/lbf_to_N)
    @printf("Wfuel   + %10.1f N (%8.1f lb)\n", parg[igWfuel], parg[igWfuel]/lbf_to_N)
    @printf("--------------------\n")
    printstyled(@sprintf("WMTO    = %10.1f N (%8.1f lb)\n\n",
                         parg[igWMTO], parg[igWMTO]/lbf_to_N); color=:bold)

    @printf("Wfuse   + %10.1f N (%8.1f lb)\n", parg[igWfuse ], parg[igWfuse ]/lbf_to_N)
    @printf("Wwing   + %10.1f N (%8.1f lb)\n", parg[igWwing ], parg[igWwing ]/lbf_to_N)
    @printf("Wvtail  + %10.1f N (%8.1f lb)\n", parg[igWvtail], parg[igWvtail]/lbf_to_N)
    @printf("Whtail  + %10.1f N (%8.1f lb)\n", parg[igWhtail], parg[igWhtail]/lbf_to_N)
    @printf("Wtesys  + %10.1f N (%8.1f lb)\n", parg[igWtesys], parg[igWtesys]/lbf_to_N)
    @printf("Wftank  + %10.1f N (%8.1f lb)\n", parg[igWftank], parg[igWftank]/lbf_to_N)
    @printf("Wadd    + %10.1f N (%8.1f lb)\n", Wtotadd, Wtotadd/lbf_to_N)
    @printf("--------------------\n")
    printstyled(@sprintf("Wempty  = %10.1f N (%8.1f lb)\n\n", 
    parg[igWfuse] + parg[igWwing]+ parg[igWvtail] + parg[igWhtail] + 
    parg[igWtesys] + +parg[igWftank] + Wtotadd, 
    (parg[igWfuse] + parg[igWwing]+ parg[igWvtail] + parg[igWhtail] + 
    parg[igWtesys] + +parg[igWftank] + Wtotadd)/lbf_to_N); color=:bold)

    @printf("Wtshaft + %10.1f N × %d\n", parg[igWtshaft], parpt[ipt_nTshaft])
    @printf("Wcat    + %10.1f N × %d\n", parg[igWcat   ], parpt[ipt_nTshaft])
    @printf("Wgen    + %10.1f N × %d\n", parg[igWgen   ], parpt[ipt_ngen]) 
    @printf("Winv    + %10.1f N × %d\n", parg[igWinv   ], parpt[ipt_nfan]) 
    @printf("Wmot    + %10.1f N × %d\n", parg[igWmot   ], parpt[ipt_nfan]) 
    @printf("Wfan    + %10.1f N × %d\n", parg[igWfan   ], parpt[ipt_nfan]) 
    @printf("--------------------\n")
    printstyled(@sprintf("Wtesys  = %10.1f N (%8.1f lb)\n\n",
     parg[igWtesys], parg[igWtesys]/lbf_to_N ), color = :bold)

    @printf("Wftnkins  = %10.1f N (%8.1f lb)\n", parg[igWinsftank ], parg[igWinsftank ]/lbf_to_N) 
    @printf("Wftank    = %10.1f N (%8.1f lb)\n\n", parg[igWftank    ], parg[igWftank    ]/lbf_to_N) 
    @printf("lftank    = %10.1f m (%8.1f ft)\n"  , parg[iglftank    ], parg[iglftank    ]/ft_to_m) 
    @printf("Rftank    = %10.1f m (%8.1f ft)\n"  , parg[igRftank    ], parg[igRftank    ]/ft_to_m) 

    @printf("ηtank   = %3.1f %% \n\n", parg[igWfuel]/(parg[igWfuel] + parg[igWftank])*100)
end

"""
`geometry` prints out the layout of the aircraft
"""
function geometry(parg)

    printstyled("Geometry Layout:\n -------------- \n", color=:bold )
    @printf("xnose   = %5.1f m (%8.1f ft)\n", parg[igxnose  ] , parg[igxnose   ]/ft_to_m)
    @printf("xend    = %5.1f m (%8.1f ft)\n", parg[igxend   ] , parg[igxend    ]/ft_to_m)
    @printf("xwing   = %5.1f m (%8.1f ft)\n", parg[igxwing  ] , parg[igxwing   ]/ft_to_m)
    @printf("xhtail  = %5.1f m (%8.1f ft)\n", parg[igxhtail ] , parg[igxhtail  ]/ft_to_m)
    @printf("xvtail  = %5.1f m (%8.1f ft)\n", parg[igxvtail ] , parg[igxvtail  ]/ft_to_m)
    @printf("xblend1 = %5.1f m (%8.1f ft)\n", parg[igxblend1] , parg[igxblend1 ]/ft_to_m)
    @printf("xblend2 = %5.1f m (%8.1f ft)\n", parg[igxblend2] , parg[igxblend2 ]/ft_to_m)
    @printf("xshell1 = %5.1f m (%8.1f ft)\n", parg[igxshell1] , parg[igxshell1 ]/ft_to_m)
    @printf("xshell2 = %5.1f m (%8.1f ft)\n", parg[igxshell2] , parg[igxshell2 ]/ft_to_m)
    @printf("xhbend  = %5.1f m (%8.1f ft)\n", parg[igxhbend ] , parg[igxhbend  ]/ft_to_m)
    @printf("xvbend  = %5.1f m (%8.1f ft)\n", parg[igxvbend ] , parg[igxvbend  ]/ft_to_m)
    @printf("xwbox   = %5.1f m (%8.1f ft)\n", parg[igxwbox  ] , parg[igxwbox   ]/ft_to_m)
    @printf("xhbox   = %5.1f m (%8.1f ft)\n", parg[igxhbox  ] , parg[igxhbox   ]/ft_to_m)
    @printf("xvbox   = %5.1f m (%8.1f ft)\n", parg[igxvbox  ] , parg[igxvbox   ]/ft_to_m)
    @printf("xtshaft = %5.1f m (%8.1f ft)\n", parg[igxtshaft] , parg[igxtshaft ]/ft_to_m)
    @printf("xgen    = %5.1f m (%8.1f ft)\n", parg[igxgen   ] , parg[igxgen    ]/ft_to_m)
    @printf("xcat    = %5.1f m (%8.1f ft)\n", parg[igxcat   ] , parg[igxcat    ]/ft_to_m)
    @printf("xftank  = %5.1f m (%8.1f ft)\n", parg[igxftank ] , parg[igxftank  ]/ft_to_m)

end