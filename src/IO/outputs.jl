"""
`weight_buildup` prints out the weight build up for the aircraft
"""
function weight_buildup(ac::aircraft; io=stdout)
    parg = ac.parg
    Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
    Whpesys = parg[igWMTO] * parg[igfhpesys]
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
    
    @printf(io, "\nRfuse  = %5.1f m (%8.1f ft)\n", parg[igRfuse ] , parg[igRfuse  ]/ft_to_m)

    
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

"""
    stickfig(parg,para,pari,parm; ax = nothing, label_fs = 16)

`stickfig` plots a "stick figure" airplane on the axis `ax` if provided 
(else a new axis is created). This is used in conjuction with other functions
like [`plot_details`](@ref) to create summary plots to track progress of optimization
or present results.
"""
function stickfig(ac::aircraft; ax = nothing, label_fs = 16)
    pari = ac.pari
    parg = ac.parg
    @views pare = ac.pare[:,:,1]
    @views para = ac.para[:,:,1]
    @views parm = ac.parm[:,:,1]
    # Wing
        co = parg[igco]
        cs = parg[igco]*parg[iglambdas]
        ct = parg[igco]*parg[iglambdat]

        sweep = parg[igsweep  ]
        λs = parg[iglambdas]
        λt = parg[iglambdat]

        bo = parg[igbo]
        bs = parg[igbs]
        b  = parg[igb ]

        xax = 0.40
        xcLE = -xax
        xcTE = 1.0 - xax

        dx = parg[igxwbox]      
        etas = bs/b
        etao = bo/b
        cosL = cos(sweep*pi/180.0)
        tanL = tan(sweep*pi/180.0)

        xs = tanL*(bs-bo)/2.0
        xt = tanL*(b -bo)/2.0

        xw = zeros(6)
        yw = zeros(6)

        #X locations of wing vertices
        xw[1] =      co*xcLE + dx
        xw[2] = xs + cs*xcLE + dx
        xw[3] = xt + ct*xcLE + dx
        xw[4] = xt + ct*xcTE + dx
        xw[5] = xs + cs*xcTE + dx
        xw[6] =      co*xcTE + dx
        
        #Y locations of wing vertices
        yw[1] = bo/2.0
        yw[2] = bs/2.0
        yw[3] = b /2.0
        yw[4] = b /2.0
        yw[5] = bs/2.0
        yw[6] = bo/2.0
    # Fuse
        Rfuse = parg[igRfuse]
        wfb   = parg[igwfb]

        anose    = parg[iganose]
        btail    = parg[igbtail]

        xnose    = parg[igxnose   ]
        xend     = parg[igxend    ]
        xblend1  = parg[igxblend1 ]
        xblend2  = parg[igxblend2 ]
        xhtail   = parg[igxhtail  ]
        xvtail   = parg[igxvtail  ]
        xwing    = parg[igxwing   ]

        xwbox    = parg[igxwbox   ]
        xhbox    = parg[igxhbox   ]
        xvbox    = parg[igxvbox   ]

        lcyl = xblend2 - xblend1
        xtail = xvtail 
        
        hwidth = Rfuse + wfb
        
        nnose = 10
        ntail = 10

        xf = zeros(nnose + ntail + 1)
        yf = zeros(nnose + ntail + 1)

        if pari[iifclose] == 0
            dytail = -hwidth 
        else
            dytail = -0.2*hwidth
        end

        for i = 1: nnose
            fraci = float(i-1)/float(nnose-1)
            fracx = cos(0.5*pi*fraci)

            k = i
            xf[k] = xblend1 + (xnose-xblend1)*fracx
            yf[k] = hwidth*(1.0 - fracx^anose)^(1.0/anose)
        end

        for i = 1: ntail
            fraci = float(i-1)/float(ntail-1)
            fracx = fraci

            k = i + nnose
            xf[k] = xblend2 + (xend-xblend2)*fracx
            yf[k] = hwidth + dytail*fracx^btail
        end

        k = k+1
        xf[k] = xf[k-1]
        yf[k] = 0.

    
    # Tail
        xh = zeros(6)
        yh = zeros(6)
        
        boh = parg[igboh]
        Sh  = parg[igSh]
        ARh = parg[igARh]
        lambdah = parg[iglambdah]
        sweeph  = parg[igsweeph]

        bh = sqrt(Sh*ARh)
        coh = Sh/(boh + (bh-boh)*0.5*(1.0+lambdah))


        dx = xhbox
        tanLh = tan(sweeph*π/180.0)
        cth = coh*lambdah

        xaxh = 0.40

        xoLEh = coh*(    - xaxh) + dx
        xoTEh = coh*(1.0 - xaxh) + dx
        xtLEh = cth*(    - xaxh) + dx + 0.5*(bh - boh)*tanLh
        xtTEh = cth*(1.0 - xaxh) + dx + 0.5*(bh - boh)*tanLh

        yoLEh = 0.5*boh
        yoTEh = 0.5*boh
        ytLEh = 0.5*bh
        ytTEh = 0.5*bh


        if (pari[iifclose] == 0)
            xcLEh = xoLEh
            xcTEh = xoTEh
            ycLEh = yoLEh
            ycTEh = yoTEh
        else
            xcLEh = coh*(    - xaxh) + dx + 0.5*(0. - boh)*tanLh
            xcTEh = coh*(1.0 - xaxh) + dx + 0.5*(0. - boh)*tanLh
            ycLEh = 0.0
            ycTEh = 0.0
        end

        if  (false)
            yoLEh = hwidth
            xoLEh = yoLEh*tanLh + xcLEh
            if (xoLEh > xblend2)
                fracx = min( (xoLEh-xblend2)/(xend-xblend2) , 1.0 )
                yoLEh = hwidth + dytail*fracx^btail
                xoLEh = yoLEh*tanLh + xcLEh
            end
            if (xoLEh > xblend2)
                fracx = min( (xoLEh-xblend2)/(xend-xblend2) , 1.0 )
                yoLEh = hwidth + dytail*fracx^btail
                xoLEh = yoLEh*tanLh + xcLEh
            end

            yoTEh = 0.0
            xoTEh = yoTEh*tanLh + xcTEh
            if (xoTEh > xblend2)
                fracx = min( (xoTEh-xblend2)/(xend-xblend2) , 1.0 )
                yoTEh = hwidth + dytail*fracx^btail
                xoTEh = yoTEh*tanLh + xcTEh
            end
            if (xoTEh > xblend2)
                fracx = min( (xoTEh-xblend2)/(xend-xblend2) , 1.0 )
                yoTEh = hwidth + dytail*fracx^btail
                xoTEh = yoTEh*tanLh + xcTEh
            end
        
        end
        


        xh[ 1] = xcLEh
        xh[ 2] = xoLEh
        xh[ 3] = xtLEh
        xh[ 4] = xtTEh
        xh[ 5] = xoTEh
        xh[ 6] = xcTEh
  
        yh[ 1] = ycLEh
        yh[ 2] = yoLEh
        yh[ 3] = ytLEh
        yh[ 4] = ytTEh
        yh[ 5] = yoTEh
        yh[ 6] = ycTEh
  

    # Fuel tank
        Rtank = Rfuse - 0.1 # Account for clearance_fuse
        l = max(parg[iglftankin], parg[iglftank])
        ARtank = 2.0
        xcyl0 = parg[igxftank] - l/2 + Rtank/ARtank
        xcyl1 = parg[igxftank] + l/2 - Rtank/ARtank
        ntank = 8
        xt = zeros(ntank*2 )
        yt = zeros(ntank*2 )
        for i = 1: ntank
            fraci = float(i-1)/float(ntank-1)
            fracx = cos(0.5*pi*fraci)

            k = i
            xt[k] = xcyl0 - Rtank/ARtank*fracx
            yt[k] = sqrt(Rtank^2 * max((1 - ((xt[k]-xcyl0)/(Rtank/ARtank))^2), 0.0) )
        end
        # k = k+1
        # xt[k] = xcyl0 + parg[iglftank]
        # yt[k] = Rtank
        for i = 1: ntank
            fraci = float(i-1)/float(ntank-1)
            fracx = sin(0.5*pi*fraci)

            k = i + ntank
            xt[k] = xcyl1 + (xcyl1 + Rtank/ARtank - xcyl1)*fracx
            yt[k] = sqrt(Rtank^2 * max((1 - ((xt[k]-xcyl1)/(Rtank/ARtank))^2), 0.0) )
        end

        # xt = LinRange(xcyl0 - Rfuse/ARtank , xcyl0, 20 )
        # yt = zeros(length(xt))
        # @. yt = sqrt(Rfuse^2 * max((1 - ((xt-xcyl0)/(Rfuse/ARtank))^2), 0.0) )

        xshell = zeros(ntank)
        yshell = zeros(ntank)
        AR = 3.0
        xshellcenter = parg[igxshell2] - Rfuse/AR
        for i = 1: ntank
            fraci = float(i-1)/float(ntank-1)
            fracx = sin(0.5*pi*fraci)

            k = i
            xshell[k] = xshellcenter + Rfuse/AR *fracx
            yshell[k] = sqrt(Rfuse^2 * max((1 - ((xshell[k]-xshellcenter)/(Rfuse/AR))^2), 0.0) )
        end

    pax = parg[igWpay]/parm[imWperpax]
    seat_pitch = 30.0 * in_to_m
    seat_width = 19.0 * in_to_m
    aisle_halfwidth = 10.0 * in_to_m # per CFR § 25.815 

    seats_per_row = Int(2*parg[igRfuse] ÷ (seat_width + aisle_halfwidth/3))
    rows = Int(ceil(pax / seats_per_row))

    println("Seats per row = $seats_per_row, Total rows = $rows")
    yseats, symmetric_seats = arrange_seats(seats_per_row, parg[igRfuse])

    xseats = zeros(rows)'
    xseats[1] = parg[igxshell1 ] + 10.0*ft_to_m 
    for r in 2:rows
        emergency_exit = 0.0
        if (r == 12 || r == 13)
            emergency_exit = seat_pitch/2
        end
        xseats[r] = xseats[r-1] + seat_pitch + emergency_exit
    end

    ## Plot
    if ax === nothing
        # plt.style.use(["../miscellaneous/prash.mplstyle"]) # HACK
        fig, ax = plt.subplots(figsize=(8,5), dpi = 300)
    else
        ax.cla()
    end
        # Plot wing
            ax.plot(xw, yw, "-k")
            ax.plot(xw, -yw, "-k")
            
            # Panel break
            ax.plot(xw[[2,5]],  yw[[2,5]], "-k", lw = 1.0, alpha = 0.5)
            ax.plot(xw[[2,5]], -yw[[2,5]], "-k", lw = 1.0, alpha = 0.5)
        
        # Plot Tail
        tailz = 1
        parg[igzhtail] > 0 ? tailz = 21 : tailz = 1
            # ax.plot(xh,  yh, "-k", zorder = tailz)
            # ax.plot(xh, -yh, "-k", zorder = tailz)
            ax.fill_between(xh, -yh, yh, facecolor = "w", alpha = 0.8, edgecolor = "k", zorder = tailz, linewidth = 2.0)
        xvt = [-0.4, -0.15, 0.2, 0.6].*parg[igcov] .+ parg[igxvbox]
        yvt = hcat([0.0 ones(length(xvt) - 2)' .*(parg[igcov]*parg[ighboxv]/2) 0.0])[:]
        ax.fill_between(xvt, -yvt, yvt, facecolor = "k", alpha = 0.8, edgecolor = "k", zorder = 22)

        # Plot fuse
            # ax.fill(xf,  yf, facecolor = "w", edgecolor = "k")
            # ax.fill(xf, -yf, facecolor = "w", edgecolor = "k")
            ax.fill_between(xf, -yf, yf, facecolor = "w", edgecolor = "k", zorder = 5, linewidth = 2.0)
            
        # Tank
        if (pari[iifwing] == 0)
            ax.plot(xt,  yt, "k", lw = 1.5, zorder = 10)
            ax.plot(xt, -yt, "k", lw = 1.5, zorder = 10)
            ax.fill_between(xt, -yt, yt, facecolor = "r", alpha = 0.1, edgecolor = "k", zorder = 6, linewidth = 1.0)
            ax.text(parg[igxftank], 0.0, "LH\$_2\$", fontsize = label_fs-2.0, zorder = 10, ha="center", va="center")
        end

        # Xshell2
        ax.plot(xshell,  yshell, "k", lw = 1.5, zorder = 10)
        ax.plot(xshell, -yshell, "k", lw = 1.5, zorder = 10)

        # Plot Engines:
            D = parg[igdaftfan]
            
            lnace = parg[iglnaceaft]
            x = parg[igxtshaft]
            ax.plot([x,x, x+lnace, x+lnace, x], [ D/8,  D/8 + D,  D/8 + D*3/4,  D/8 + 1/4*D,  D/8], lw = 1.5, color = "r", zorder = 25)
            ax.plot([x,x, x+lnace, x+lnace, x], [-D/8, -D/8 - D, -D/8 - D*3/4, -D/8 - 1/4*D, -D/8], lw = 1.5, color = "r", zorder = 25)

            D = parg[igdfan]
            neng = parg[igneng]
            lnace = parg[iglnace]
            ηs = parg[igetas]
            dy = 2*D # space to leave near wing root and tip [m]
            if parg[igneng] == 2
                yi = [ηs*b/2]
            else
                yi = LinRange(bo/2 + dy , b/2 *3/4, Int(parg[igneng]/2))
            end
            xi = zero(yi)
            ηi = yi/(b/2)
            ηs = bs/b
            ηo = bo/b
            ci = zero(yi)
            for (i, η)  in enumerate(ηi)
                if η <=ηs
                    ci[i] = co*(1  + (λs -  1)*(η - ηo)/(ηs - ηo))
                else
                    ci[i] = co*(λs + (λt - λs)*(η - ηs)/(1  - ηs))
                end
            end

            tanL = tan(parg[igsweep]*π/180.0)
            @. xi = tanL * (yi - bo/2) - 0.4ci + parg[igxwbox] - 1.0
            ax.plot( [xi, xi, xi.+lnace, xi.+lnace, xi] , [yi.-D/2, yi.+D/2, yi.+D/3, yi.-D/3, yi.-D/2 ], color = "r", lw = 1.5)

        # Plot NP and CG range
            ax.scatter(parg[igxNP], 0.0, color = "k", marker="o", zorder = 21, label = "NP")
            ax.text(parg[igxNP], -1.0, "NP", fontsize=label_fs-2.0, ha="center", va="center", zorder = 21)

            ax.annotate("", xy=(parg[igxCGfwd ] , 0.0), xytext=(parg[igxCGaft ] , 0.0),
            fontsize=16, ha="center", va="bottom",
            arrowprops=Dict("arrowstyle"=> "|-|, widthA=0.2, widthB=0.2"),
            zorder = 21, label = "CG movement")
            ax.text(0.5*(parg[igxCGfwd ]+parg[igxCGaft ]), -1.0, "CG", fontsize=label_fs-2.0, ha="center", va="center", zorder = 21)

        # Show seats
        if symmetric_seats
            ax.scatter(ones(length(yseats),1).*xseats, ones(1,rows).* yseats, color = "gray", alpha = 0.1, marker = "s", s=15, zorder = 21)
            ax.scatter(ones(length(yseats),1).*xseats, ones(1,rows).*-yseats, color = "gray", alpha = 0.1, marker = "s", s=15, zorder = 21)
        else
            ax.scatter(ones(length(yseats),1).*xseats, ones(1,rows).* yseats, color = "gray", alpha = 0.1, marker = "s", s=15, zorder = 21)
        end
     # diagnostic marks
    #  ax.scatter(parg[igxftank] - l/2, 0.0, color = "k", marker="o", zorder = 21)
    #  ax.scatter(parg[igxftank], 0.0, color = "b", marker="o", zorder = 21)
    #  ax.scatter(parg[igxblend2], 0.0, color = "k", marker="o", zorder = 21)
    #  ax.plot([parg[igxftank]-l/2, parg[igxftank]+l/2],[0.0, 0.0], zorder = 21)



    # Annotations
    ax.text(0, 16, @sprintf("PFEI = %5.3f J/Nm\nM\$_{cruise}\$ = %.2f\nWMTO = %.1f tonnes\nSpan = %5.1f m\nco    = %5.1f m\n\$ \\Lambda \$ = %.1f\$^\\circ\$\nRfuse = %5.1f m\nL/D = %3.2f",
     parm[imPFEI], para[iaMach, ipcruise1],parg[igWMTO]/9.81/1000, parg[igb], parg[igco], parg[igsweep], parg[igRfuse], para[iaCL, ipcruise1]/para[iaCD, ipcruise1]),
     fontsize = label_fs, ha="left", va="top")

    yloc = -20
    ax.annotate("", xy=(0.0, yloc), xytext=( xf[end], yloc),
            fontsize=16, ha="center", va="bottom",
            arrowprops=Dict("arrowstyle"=> "|-|, widthA=0.5, widthB=0.5"),
             zorder = 30)
    ax.text(xend/2, yloc, @sprintf("l = %5.1f m", xend), bbox=Dict("ec"=>"w", "fc"=>"w"), ha="center", va="center", fontsize = 14, zorder = 31)

    # Span annotations:
    codeD = true
    codeE = false
    xcodeD = -2.0
    xcodeE = -3.5
        if codeD
            # ICAO code D 
            bmaxD = 36
            ax.vlines(xcodeD, -bmaxD/2, bmaxD/2, lw = 5, alpha = 0.2, color = "y")
            ax.hlines( bmaxD/2, xcodeD, 40.0, lw = 5, alpha = 0.2, color = "y")
            ax.hlines(-bmaxD/2, xcodeD, 40.0, lw = 5, alpha = 0.2, color = "y")
            ax.text(20, bmaxD/2+1, "ICAO Code D/ FAA Group III", color = "y", alpha = 0.8, fontsize = 12, ha="center", va="center")
        end
        if codeE
            # ICAO code E
            bmaxE = 52
            ax.vlines(xcodeE, -bmaxE/2, bmaxE/2, lw = 5, alpha = 0.2, color = "b")
            ax.hlines( bmaxE/2, xcodeE, 40.0, lw = 5, alpha = 0.2, color = "b")
            ax.hlines(-bmaxE/2, xcodeE, 40.0, lw = 5, alpha = 0.2, color = "b")
            ax.text(20, bmaxE/2+1, "ICAO Code E/ FAA Group IV", color = "b", alpha = 0.5, fontsize = 12, ha="center", va="center")
        end

    if codeE
        ax.set_ylim(min(-27, -b/2.0), max(27, b/2.0))
    elseif codeD
        ax.set_ylim(min(-23, -b/2.0),max(23, b/2.0))
    else
        ax.set_ylim(min(-20, -2b/2.0), max(20, b/2.0))
    end
    ax.set_aspect(1)
    ax.set_ylabel("y[m]")
    ax.set_xlabel("x[m]")
    plt.tight_layout()
    # ax.legend()

    return ax
end

"""
    plot_details(parg, pari, para, parm; ax = nothing)

`plot_details` combines a [`stickfig`](@ref) plot along with a mission summary,
weight and drag buildup stacked bar charts to present results.
"""
function plot_details(ac::aircraft; ax = nothing)

    pari = ac.pari
    parg = ac.parg
    @views pare = ac.pare[:,:,1]
    @views para = ac.para[:,:,1]
    @views parm = ac.parm[:,:,1]
        ## Create empty plot
        if ax === nothing
            # plt.style.use(["../miscellaneous/prash.mplstyle", "seaborn-colorblind"]) # HACK
            fig, atemp = plt.subplots(2, 2, figsize=(8,5), dpi = 300, gridspec_kw=Dict("height_ratios"=>[1, 3], "width_ratios"=>[1,3]))
            gs = atemp[1,2].get_gridspec()
            gssub = matplotlib.gridspec.SubplotSpec(gs, 0,1)
            atemp[1,1].remove()
            atemp[1,2].remove()
            axbig = fig.add_subplot(gssub)
            ax = [axbig, atemp[2,1], atemp[2,2]]
        else
            for a in ax
                a.cla()
            end
        end

        # Drag build-up
        LoD     = para[iaCL, ipcruise1]/ para[iaCD, ipcruise1]
        CL      = para[iaCL, ipcruise1]
        CD      = para[iaCD, ipcruise1]
        CDfuse  = para[iaCDfuse, ipcruise1]
        CDi     = para[iaCDi, ipcruise1]
        CDwing  = para[iaCDwing, ipcruise1]
        CDhtail = para[iaCDhtail, ipcruise1]
        CDvtail = para[iaCDvtail, ipcruise1]
        CDnace  = para[iaCDnace, ipcruise1]

        CDfusefrac  = CDfuse /CD
        CDifrac     = CDi    /CD
        CDwingfrac  = CDwing /CD
        CDhtailfrac = CDhtail/CD
        CDvtailfrac = CDvtail/CD
        CDnacefrac  = CDnace /CD

        # Weight build-up
        Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
        Whpesys = parg[igWMTO] * parg[igfhpesys]
        Wlgnose = parg[igWMTO] * parg[igflgnose]
        Wlgmain = parg[igWMTO] * parg[igflgmain]
        Wtotadd = Whpesys + Wlgnose + Wlgmain
        
        Wpay  = parg[igWpay]
        Wfuel = parg[igWfuel]
        WMTO  = parg[igWMTO]

        Wwing  = parg[igWwing]
        Wfuse  = parg[igWfuse]
        Wvtail = parg[igWvtail]
        Whtail = parg[igWhtail]
        Wtesys = parg[igWtesys]
        Wftank = parg[igWftank]

        Wwingfrac = Wwing /WMTO
        Wfusefrac = Wfuse /WMTO
        Wvtailfrac = Wvtail/WMTO
        Whtailfrac = Whtail/WMTO
        Wtesysfrac = Wtesys/WMTO
        Wftankfrac = Wftank/WMTO
        Wtotaddfrac = Wtotadd/WMTO

        Wemptyfrac = Wempty/WMTO
        Wfuelfrac  = Wfuel /WMTO
        Wpayfrac   = Wpay  /WMTO

        a = ax[2]
        bar_width = 0.2
        # ax[1,1].bar("CL", CL)
        CDbars = []
        push!(CDbars, a.bar(0, CDifrac    , width = bar_width, bottom = CDfusefrac+CDwingfrac+CDhtailfrac+CDvtailfrac+CDnacefrac, label = "CDi"))
        push!(CDbars, a.bar(0, CDnacefrac , width = bar_width, bottom = CDfusefrac+CDwingfrac+CDhtailfrac+CDvtailfrac           , label = "CDnace"))
        push!(CDbars, a.bar(0, CDvtailfrac, width = bar_width, bottom = CDfusefrac+CDwingfrac+CDhtailfrac                       , label = "CDvtail"))
        push!(CDbars, a.bar(0, CDhtailfrac, width = bar_width, bottom = CDfusefrac+CDwingfrac                                   , label = "CDhtail"))
        push!(CDbars, a.bar(0, CDwingfrac , width = bar_width, bottom = CDfusefrac                                              , label = "CDwing"))
        push!(CDbars, a.bar(0, CDfusefrac , width = bar_width, label = "CDfuse"))
        
        CDlabels = ["CDi", "CDnace", "CDvtail", "CDhtail", "CDwing", "CDfuse"]

        label_bars(a, CDbars, CDlabels; val_multiplier = CD)
        # a.legend(loc = "upper center")
        # a.legend(bbox_to_anchor=(1.05, 1))
        a.set_xlim(-1,3.5)

        
        Wbar1 = a.bar(1.5, Wpayfrac  , bottom = Wemptyfrac + Wfuelfrac, width = bar_width, label = "Wpay")
        Wbar2 = a.bar(1.5, Wfuelfrac , bottom = Wemptyfrac, width = bar_width, label = "Wfuel")
        Wbar3 = a.bar(1.5, Wemptyfrac, width = bar_width, label = "Wempty")
        Wbars = [Wbar1, Wbar2, Wbar3]
        Wlabels = ["Wpay", "Wfuel", "Wempty"]
        label_bars(a, Wbars, Wlabels, val_multiplier = WMTO/9.81/1000)
        
        Webars = []
        push!(Webars, a.bar(3, Wtotaddfrac , bottom = Wfusefrac+Wwingfrac+Whtailfrac+Wvtailfrac+Wtesysfrac+Wftankfrac, width = bar_width, label = "Wadd"))
        push!(Webars, a.bar(3, Wftankfrac  , bottom = Wfusefrac+Wwingfrac+Whtailfrac+Wvtailfrac+Wtesysfrac, width = bar_width, label = "Wftank"))
        push!(Webars, a.bar(3, Wtesysfrac  , bottom = Wfusefrac+Wwingfrac+Whtailfrac+Wvtailfrac, width = bar_width, label = "Wtesys"))
        push!(Webars, a.bar(3, Wvtailfrac  , bottom = Wfusefrac+Wwingfrac+Whtailfrac, width = bar_width, label = "Wvtail"))
        push!(Webars, a.bar(3, Whtailfrac  , bottom = Wfusefrac+Wwingfrac, width = bar_width, label = "Whtail"))
        push!(Webars, a.bar(3, Wwingfrac   , bottom = Wfusefrac, width = bar_width, label = "Wwing"))
        push!(Webars, a.bar(3, Wfusefrac   , width = bar_width, label = "Wfuse"))
        
        Welabels = ["Wadd" "Wftank" "Wtesys" "Wvtail" "Whtail" "Wwing" "Wfuse"]
        label_bars(a, Webars, Welabels, val_multiplier = WMTO/9.81/1000)

        a.hlines(Wemptyfrac, 1.5+bar_width/2, 3-bar_width/2, lw=0.8, color = "k", ls = "--")
        a.grid()
        a.set_xticks([0, 1.5, 3])
        a.set_xticklabels(["CD","WMTO", "Wempty"])
        # ar.cla()
        # ar = a.twinx()
        # ar.bar(1.7, CDi     , width = 0.4, bottom = CDfuse +CDwing +CDhtail +CDvtail +CDnace , label = "CDi")
        # ar.bar(1.7, CDnace  , width = 0.4, bottom = CDfuse +CDwing +CDhtail +CDvtail            , label = "CDnace")
        # ar.bar(1.7, CDvtail , width = 0.4, bottom = CDfuse +CDwing +CDhtail                        , label = "CDvtail")
        # ar.bar(1.7, CDhtail , width = 0.4, bottom = CDfuse +CDwing                                    , label = "CDhtail")
        # ar.bar(1.7, CDwing  , width = 0.4, bottom = CDfuse                                               , label = "CDwing")
        # ar.bar(1.7, CDfuse  , width = 0.4, label = "CDfuse")
        # ar.grid()



        # Draw mission profile
        a = ax[1]
        h     = [para[iaalt,ipclimb1:ipcruisen]; 0.0]./ft_to_m./1000 # show in 1000s of ft.
        R     = [para[iaRange,ipclimb1:ipcruisen]; para[iaRange, ipdescentn]]./nmi_to_m
        deNOx = pare[iedeNOx, :]
        fracW = [para[iafracW, ipclimb1:ipcruisen]; para[iafracW, ipdescentn]]
        mdotf = pare[iemdotf, :]
        mdotH2O = pare[iemdotf, :].*9.0
        gamV = [para[iagamV, ipclimb1:ipcruisen]; para[iagamV, ipdescentn]]

        a.plot(R, h)
        a.set_ylim(0, 60.0)
        a.set_xlabel("Range [nmi]")
        a.set_ylabel("Altitude [kft]")
        
        # ar = a.twinx()
        # ar.plot(R, gamV, color = "r")
        # ar.axhline(0.015, lw = 1.0, color = "r")
        # ar.set_ylabel("Climb angle")

        # Draw stick figure to keep track
        stickfig(ac; ax = ax[3], label_fs = 12)
        plt.tight_layout()

        #Print other details:
        # ax[3].text(48,20, @sprintf("WMTO = %.1f tons\n\$ \\Lambda \$ = %.1f\$ ^\\circ \$\n", parg[igWMTO]/9.81/1000, parg[igsweep]), va = "top")

        return ax

end

"""
Simple utility function to label bars in a stacked bar chart
"""
function label_bars(a, Bararray, labels; val_multiplier = 1, fontsize = 8)
    for (i,bar) in enumerate(Bararray)
        w, h = bar[1].get_width(), bar[1].get_height()
        x, y = bar[1].get_x(), bar[1].get_y()
        a.text(x+w, y+h/2, @sprintf("%7.3f", h*val_multiplier), ha = "left", va = "center", fontsize = fontsize)
        a.text(x-w/2, y+h/2, @sprintf("%7s", labels[i]), ha = "right", va = "center", fontsize = fontsize)
    end
end

#737-800 from TASOPT w 220 pax
function plot737compare(;weightdetail= true, fracs = false)
    fig, ax = plt.subplots(1,2,figsize=(8,5), dpi = 300)
    
    Wempty  = 105757.5* lbf_to_N

    Wpay    =  47318.9* lbf_to_N
    Wfuel   =  49862.7* lbf_to_N
    WMTO    = 202939.1* lbf_to_N

    Wtotadd =  13191.0* lbf_to_N
    Wfuse   =  44496.3* lbf_to_N
    Wwing   =  26043.3* lbf_to_N
    Whtail  =   2650.0* lbf_to_N
    Wvtail  =   1749.0* lbf_to_N
    Weng    =  17628.0* lbf_to_N

    if fracs
        Wemptyfrac  = Wempty  /WMTO
        Wpayfrac    = Wpay    /WMTO
        Wfuelfrac   = Wfuel   /WMTO
        Wtotaddfrac = Wtotadd /WMTO
        Wfusefrac   = Wfuse   /WMTO
        Wwingfrac   = Wwing   /WMTO
        Whtailfrac  = Whtail  /WMTO
        Wvtailfrac  = Wvtail  /WMTO
        Wengfrac    = Weng    /WMTO
    else
        Wemptyfrac  = Wempty  /9.81/1000 
        Wpayfrac    = Wpay    /9.81/1000 
        Wfuelfrac   = Wfuel   /9.81/1000 
        Wtotaddfrac = Wtotadd /9.81/1000
        Wfusefrac   = Wfuse   /9.81/1000
        Wwingfrac   = Wwing   /9.81/1000
        Whtailfrac  = Whtail  /9.81/1000
        Wvtailfrac  = Wvtail  /9.81/1000
        Wengfrac    = Weng    /9.81/1000
    end

    bar_width = 0.2
    a = ax[1]
    if weightdetail == false
        Wbar1 = a.bar(0.0, Wpayfrac  , color = "#0072B2", bottom = Wemptyfrac + Wfuelfrac, width = bar_width, label = "Wpay")
        Wbar2 = a.bar(0.0, Wfuelfrac , color = "#009E73", bottom = Wemptyfrac, width = bar_width, label = "Wfuel")
        Wbar3 = a.bar(0.0, Wemptyfrac, color = "#D55E00", width = bar_width, label = "Wempty")
        Wbars = [Wbar1, Wbar2, Wbar3]
        Wlabels = ["", "", ""]
        if fracs
            label_bars(a, Wbars, Wlabels, val_multiplier = WMTO/9.81/1000, fontsize = 18)
        else
            label_bars(a, Wbars, Wlabels, val_multiplier = 1/(WMTO/9.81/1000), fontsize = 18)
        end
        # =======
        #   ZIA
        # =======
        Wpay  = parg[igWpay]
        Wfuel = parg[igWfuel]
        WMTO  = parg[igWMTO]
        Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
        if fracs
            Wemptyfrac = Wempty/WMTO
            Wfuelfrac  = Wfuel /WMTO
            Wpayfrac   = Wpay  /WMTO
        else
            Wemptyfrac = Wempty/9.81/1000 #/WMTO
            Wfuelfrac  = Wfuel /9.81/1000 #/WMTO
            Wpayfrac   = Wpay  /9.81/1000 #/WMTO
        end    

        Wbar1 = a.bar(1.0, Wpayfrac  ,  color = "#0072B2", bottom = Wemptyfrac + Wfuelfrac, width = bar_width, label = "Wpay")
        Wbar2 = a.bar(1.0, Wfuelfrac ,  color = "#009E73", bottom = Wemptyfrac, width = bar_width, label = "Wfuel")
        Wbar3 = a.bar(1.0, Wemptyfrac,  color = "#D55E00", width = bar_width, label = "Wempty")
        Wbars = [Wbar1, Wbar2, Wbar3]
        Wlabels = ["Wpay", "Wfuel", "Wempty"]
        if fracs
            label_bars(a, Wbars, Wlabels, val_multiplier = WMTO/9.81/1000, fontsize = 18)
        else
            label_bars(a, Wbars, Wlabels, val_multiplier = 1/(WMTO/9.81/1000), fontsize = 18)
        end
    else
        
        Wbars = []
        push!(Wbars, a.bar(0.0, Wpayfrac  , color = "#0072B2", bottom = Wemptyfrac + Wfuelfrac, width = bar_width, label = "Wpay"))
        push!(Wbars, a.bar(0.0, Wfuelfrac , color = "#009E73", bottom = Wemptyfrac, width = bar_width, label = "Wfuel"))
        push!(Wbars, a.bar(0.0, Wtotaddfrac, width = bar_width, bottom = Wfusefrac + Wwingfrac + Whtailfrac + Wvtailfrac + Wengfrac, label = "Wadd"))
        push!(Wbars, a.bar(0.0, Wengfrac, width = bar_width, bottom = Wfusefrac + Wwingfrac + Whtailfrac + Wvtailfrac, label = "Weng"))
        push!(Wbars, a.bar(0.0, Wvtailfrac, width = bar_width, bottom = Wfusefrac + Wwingfrac + Whtailfrac, label = "Wvtail"))
        push!(Wbars, a.bar(0.0, Whtailfrac, width = bar_width, bottom = Wfusefrac + Wwingfrac, label = "Whtail"))
        push!(Wbars, a.bar(0.0, Wwingfrac, width = bar_width, bottom = Wfusefrac, label = "Wwing"))
        push!(Wbars, a.bar(0.0, Wfusefrac,color = "k",  width = bar_width, label = "Wfuse"))
        Wlabels = fill("", length(Wbars))
        if fracs
            label_bars(a, Wbars, Wlabels, val_multiplier = WMTO/9.81/1000, fontsize = 18)
        else
            label_bars(a, Wbars, Wlabels, val_multiplier = 1/(WMTO/9.81/1000), fontsize = 18)
        end

        # =======
        #   ZIA
        # =======
        Wpay  = parg[igWpay]
        Wfuel = parg[igWfuel]
        WMTO  = parg[igWMTO]
        Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
        Wwing = parg[igWwing]
        Wfuse = parg[igWfuse]
        Whtail = parg[igWhtail]
        Wvtail = parg[igWvtail]
        Weng = parg[igWtesys]
            Whpesys = parg[igWMTO] * parg[igfhpesys]
            Wlgnose = parg[igWMTO] * parg[igflgnose]
            Wlgmain = parg[igWMTO] * parg[igflgmain]
        Wtotadd = Whpesys + Wlgnose + Wlgmain
        Wftank = parg[igWftank]

        if fracs
            Wemptyfrac = Wempty /WMTO
            Wfuelfrac  = Wfuel  /WMTO
            Wpayfrac   = Wpay   /WMTO
            Wtotaddfrac = Wtotadd /WMTO
            Wfusefrac   = Wfuse   /WMTO
            Wwingfrac   = Wwing   /WMTO
            Whtailfrac  = Whtail  /WMTO
            Wvtailfrac  = Wvtail  /WMTO
            Wengfrac    = Weng    /WMTO
            Wftankfrac  = Wftank  /WMTO
        else
            Wemptyfrac = Wempty/9.81/1000 #/WMTO
            Wfuelfrac  = Wfuel /9.81/1000 #/WMTO
            Wpayfrac   = Wpay  /9.81/1000 #/WMTO
            Wtotaddfrac = Wtotadd /9.81/1000
            Wfusefrac   = Wfuse   /9.81/1000
            Wwingfrac   = Wwing   /9.81/1000
            Whtailfrac  = Whtail  /9.81/1000
            Wvtailfrac  = Wvtail  /9.81/1000
            Wengfrac    = Weng    /9.81/1000
            Wftankfrac  = Wftank  /9.81/1000
        end

        Wbars = []
        push!(Wbars, a.bar(1.0, Wpayfrac  , color = "#0072B2", bottom = Wemptyfrac + Wfuelfrac, width = bar_width, label = "Wpay"))
        push!(Wbars, a.bar(1.0, Wfuelfrac , color = "#009E73", bottom = Wemptyfrac, width = bar_width, label = "Wfuel"))
        push!(Wbars, a.bar(1.0, Wftankfrac , bottom = Wfusefrac + Wwingfrac + Whtailfrac + Wvtailfrac + Wengfrac+Wtotaddfrac, width = bar_width, label = "Wftank"))
        push!(Wbars, a.bar(1.0, Wtotaddfrac, width = bar_width, bottom = Wfusefrac + Wwingfrac + Whtailfrac + Wvtailfrac + Wengfrac, label = "Wadd"))
        push!(Wbars, a.bar(1.0, Wengfrac, width = bar_width, bottom = Wfusefrac + Wwingfrac + Whtailfrac + Wvtailfrac, label = "Weng"))
        push!(Wbars, a.bar(1.0, Wvtailfrac, width = bar_width, bottom = Wfusefrac + Wwingfrac + Whtailfrac, label = "Wvtail"))
        push!(Wbars, a.bar(1.0, Whtailfrac, width = bar_width, bottom = Wfusefrac + Wwingfrac, label = "Whtail"))
        push!(Wbars, a.bar(1.0, Wwingfrac, width = bar_width, bottom = Wfusefrac, label = "Wwing"))
        push!(Wbars, a.bar(1.0, Wfusefrac, color = "k", width = bar_width, label = "Wfuse"))
        Wlabels = ["Wpay", "Wfuel", "Wftank", "Wadd", "Wprop", "Wvtail", "Whtail", "Wwing", "Wfuse" ]
        if fracs
            label_bars(a, Wbars, Wlabels, val_multiplier = WMTO/9.81/1000, fontsize = 18)
        else
            label_bars(a, Wbars, Wlabels, val_multiplier = 1/(WMTO/9.81/1000), fontsize = 18)
        end
    end
    a.set_ylabel("Mass [tonnes]", fontsize = 20)
    a.set_xticks([0, 1])
    a.set_xticklabels(["737-800", "ZIA"], fontsize = 20)
    a.grid()
    
    CD     = 0.03185
    CDi    = 0.01114
    CDfuse = 0.00622
    CDwing = 0.00830
    CDhtail= 0.00261
    CDvtail= 0.00176
    CDnace = 0.00182
    
    CDifrac    = CDi     #/CD
    CDfusefrac = CDfuse  #/CD
    CDwingfrac = CDwing  #/CD
    CDhtailfrac= CDhtail #/CD
    CDvtailfrac= CDvtail #/CD
    CDnacefrac = CDnace  #/CD
    
    a = ax[2]
    CDbars = []
    push!(CDbars, a.bar(0, CDifrac    , width = bar_width, bottom = CDfusefrac+CDwingfrac+CDhtailfrac+CDvtailfrac+CDnacefrac, label = "CDi"))
    push!(CDbars, a.bar(0, CDnacefrac , width = bar_width, bottom = CDfusefrac+CDwingfrac+CDhtailfrac+CDvtailfrac           , label = "CDnace"))
    push!(CDbars, a.bar(0, CDvtailfrac, width = bar_width, bottom = CDfusefrac+CDwingfrac+CDhtailfrac                       , label = "CDvtail"))
    push!(CDbars, a.bar(0, CDhtailfrac, width = bar_width, bottom = CDfusefrac+CDwingfrac                                   , label = "CDhtail"))
    push!(CDbars, a.bar(0, CDwingfrac , width = bar_width, bottom = CDfusefrac                                              , label = "CDwing"))
    push!(CDbars, a.bar(0, CDfusefrac , width = bar_width, label = "CDfuse"))
    
    CDlabels = ["", "", "", "", "", ""]
    
    label_bars(a, CDbars, CDlabels; val_multiplier = 1/CD, fontsize = 18)
    
    
    CD      = para[iaCD, ipcruise1]
    CDfuse  = para[iaCDfuse, ipcruise1]
    CDi     = para[iaCDi, ipcruise1]
    CDwing  = para[iaCDwing, ipcruise1]
    CDhtail = para[iaCDhtail, ipcruise1]
    CDvtail = para[iaCDvtail, ipcruise1]
    CDnace  = para[iaCDnace, ipcruise1]
    
    CDfusefrac  = CDfuse  #/CD
    CDifrac     = CDi     #/CD
    CDwingfrac  = CDwing  #/CD
    CDhtailfrac = CDhtail #/CD
    CDvtailfrac = CDvtail #/CD
    CDnacefrac  = CDnace  #/CD
    CDbars = []
    push!(CDbars, a.bar(1, CDifrac    , width = bar_width, bottom = CDfusefrac+CDwingfrac+CDhtailfrac+CDvtailfrac+CDnacefrac, label = "CDi"))
    push!(CDbars, a.bar(1, CDnacefrac , width = bar_width, bottom = CDfusefrac+CDwingfrac+CDhtailfrac+CDvtailfrac           , label = "CDnace"))
    push!(CDbars, a.bar(1, CDvtailfrac, width = bar_width, bottom = CDfusefrac+CDwingfrac+CDhtailfrac                       , label = "CDvtail"))
    push!(CDbars, a.bar(1, CDhtailfrac, width = bar_width, bottom = CDfusefrac+CDwingfrac                                   , label = "CDhtail"))
    push!(CDbars, a.bar(1, CDwingfrac , width = bar_width, bottom = CDfusefrac                                              , label = "CDwing"))
    push!(CDbars, a.bar(1, CDfusefrac , width = bar_width, label = "CDfuse"))
    
    CDlabels = ["CDi", "CDnace", "CDvtail", "CDhtail", "CDwing", "CDfuse"]
    
    label_bars(a, CDbars, CDlabels; val_multiplier = 1/CD, fontsize = 18)
    
    a.set_ylabel("CD [-]", fontsize = 20)
    a.set_xticks([0, 1])
    a.set_xticklabels(["737-800", "ZIA"], fontsize = 20)
    a.grid()
    plt.tight_layout()
    end

"""
    MomentShear(parg)

Plot moment and shear diagrams
"""
function MomentShear(parg)
    co = parg[igco]
    cs = parg[igco]*parg[iglambdas]
    ct = parg[igco]*parg[iglambdat]
  
    bo = parg[igbo]
    bs = parg[igbs]
    b  = parg[igb ]

    etas = bs/b
    etao = bo/b

    λs = parg[iglambdas]
    λt = parg[iglambdat]

    Ss = parg[igSsmax]
    Ms = parg[igMsmax]
    So = parg[igSomax]
    Mo = parg[igMomax]
    
    etaRange =[LinRange(etao,etas,20) ; LinRange(etas,1,20)]
    c = zero(etaRange)
    S = zero(etaRange)
    M = zero(etaRange)

    for (i,eta) in enumerate(etaRange)
        if eta<etao
            c[i] = co
        elseif eta<etas
            c[i] = co*( 1 + (λs - 1)*(eta - etao)/(etas - etao))
            S[i] = So
            M[i] = Mo
        else
            c[i] = co*(λs + (λt -λs)*(eta - etas)/(   1 - etas))
            S[i] = Ss*(c[i]/cs)^2
            M[i] = Ms*(c[i]/cs)^3
        end
    end
    fig, ax = plt.subplots(2,1,figsize=(8,5), sharex = true, dpi = 300)
    ax[1].plot(etaRange,S)
    ax[1].set_ylabel("Shear")
    ax[2].plot(etaRange,M)
    ax[2].set_ylabel("Moment")

    for a in ax
        a.axvline(etas, ls = "--")
        a.axvline(etao)
    end

    return ax

end

"""
    high_res_airplane_plot(parg, pari, parm; ax = nothing, label_fs = 16, save_name = nothing)

plots high resolution figure for publications
"""
function high_res_airplane_plot(ac; ax = nothing, label_fs = 16, save_name = nothing)

    pari = ac.pari
    parg = ac.parg
    @views pare = ac.pare[:,:,1]
    @views para = ac.para[:,:,1]
    @views parm = ac.parm[:,:,1]
    # Wing
        co = parg[igco]
        cs = parg[igco]*parg[iglambdas]
        ct = parg[igco]*parg[iglambdat]

        sweep = parg[igsweep  ]
        λs = parg[iglambdas]
        λt = parg[iglambdat]

        bo = parg[igbo]
        bs = parg[igbs]
        b  = parg[igb ]

        xax = 0.40
        xcLE = -xax
        xcTE = 1.0 - xax

        dx = parg[igxwbox]      
        etas = bs/b
        etao = bo/b
        cosL = cos(sweep*pi/180.0)
        tanL = tan(sweep*pi/180.0)

        xs = tanL*(bs-bo)/2.0
        xt = tanL*(b -bo)/2.0

        xw = zeros(6)
        yw = zeros(6)

        #X locations of wing vertices
        xw[1] =      co*xcLE + dx
        xw[2] = xs + cs*xcLE + dx
        xw[3] = xt + ct*xcLE + dx
        xw[4] = xt + ct*xcTE + dx
        xw[5] = xs + cs*xcTE + dx
        xw[6] =      co*xcTE + dx
        
        #Y locations of wing vertices
        yw[1] = bo/2.0
        yw[2] = bs/2.0
        yw[3] = b /2.0
        yw[4] = b /2.0
        yw[5] = bs/2.0
        yw[6] = bo/2.0
    # Fuse
        Rfuse = parg[igRfuse]
        wfb   = parg[igwfb]

        anose    = parg[iganose]
        btail    = parg[igbtail]

        xnose    = parg[igxnose   ]
        xend     = parg[igxend    ]
        xblend1  = parg[igxblend1 ]
        xblend2  = parg[igxblend2 ]
        xhtail   = parg[igxhtail  ]
        xvtail   = parg[igxvtail  ]
        xwing    = parg[igxwing   ]

        xwbox    = parg[igxwbox   ]
        xhbox    = parg[igxhbox   ]
        xvbox    = parg[igxvbox   ]

        lcyl = xblend2 - xblend1
        xtail = xvtail 
        
        hwidth = Rfuse + wfb
        
        nnose = 15
        ntail = 10

        xf = zeros(nnose + ntail + 1)
        yf = zeros(nnose + ntail + 1)

        if pari[iifclose] == 0
            dytail = -hwidth 
        else
            dytail = -0.2*hwidth
        end

        for i = 1: nnose
            fraci = float(i-1)/float(nnose-1)
            fracx = cos(0.5*pi*fraci)

            k = i
            xf[k] = xblend1 + (xnose-xblend1)*fracx
            yf[k] = hwidth*(1.0 - fracx^anose)^(1.0/anose)
        end

        for i = 1: ntail
            fraci = float(i-1)/float(ntail-1)
            fracx = fraci

            k = i + nnose
            xf[k] = xblend2 + (xend-xblend2)*fracx
            yf[k] = hwidth + dytail*fracx^btail
        end

        k = k+1
        xf[k] = xf[k-1]
        yf[k] = 0.

    
    # Tail
        xh = zeros(6)
        yh = zeros(6)
        
        boh = parg[igboh]
        Sh  = parg[igSh]
        ARh = parg[igARh]
        lambdah = parg[iglambdah]
        sweeph  = parg[igsweeph]

        bh = sqrt(Sh*ARh)
        coh = Sh/(boh + (bh-boh)*0.5*(1.0+lambdah))


        dx = xhbox
        tanLh = tan(sweeph*π/180.0)
        cth = coh*lambdah

        xaxh = 0.40

        xoLEh = coh*(    - xaxh) + dx
        xoTEh = coh*(1.0 - xaxh) + dx
        xtLEh = cth*(    - xaxh) + dx + 0.5*(bh - boh)*tanLh
        xtTEh = cth*(1.0 - xaxh) + dx + 0.5*(bh - boh)*tanLh

        yoLEh = 0.5*boh
        yoTEh = 0.5*boh
        ytLEh = 0.5*bh
        ytTEh = 0.5*bh


        if (pari[iifclose] == 0)
            xcLEh = xoLEh
            xcTEh = xoTEh
            ycLEh = yoLEh
            ycTEh = yoTEh
        else
            xcLEh = coh*(    - xaxh) + dx + 0.5*(0. - boh)*tanLh
            xcTEh = coh*(1.0 - xaxh) + dx + 0.5*(0. - boh)*tanLh
            ycLEh = 0.0
            ycTEh = 0.0
        end

        if  (false)
            yoLEh = hwidth
            xoLEh = yoLEh*tanLh + xcLEh
            if (xoLEh > xblend2)
                fracx = min( (xoLEh-xblend2)/(xend-xblend2) , 1.0 )
                yoLEh = hwidth + dytail*fracx^btail
                xoLEh = yoLEh*tanLh + xcLEh
            end
            if (xoLEh > xblend2)
                fracx = min( (xoLEh-xblend2)/(xend-xblend2) , 1.0 )
                yoLEh = hwidth + dytail*fracx^btail
                xoLEh = yoLEh*tanLh + xcLEh
            end

            yoTEh = 0.0
            xoTEh = yoTEh*tanLh + xcTEh
            if (xoTEh > xblend2)
                fracx = min( (xoTEh-xblend2)/(xend-xblend2) , 1.0 )
                yoTEh = hwidth + dytail*fracx^btail
                xoTEh = yoTEh*tanLh + xcTEh
            end
            if (xoTEh > xblend2)
                fracx = min( (xoTEh-xblend2)/(xend-xblend2) , 1.0 )
                yoTEh = hwidth + dytail*fracx^btail
                xoTEh = yoTEh*tanLh + xcTEh
            end
        
        end
        


        xh[ 1] = xcLEh
        xh[ 2] = xoLEh
        xh[ 3] = xtLEh
        xh[ 4] = xtTEh
        xh[ 5] = xoTEh
        xh[ 6] = xcTEh
  
        yh[ 1] = ycLEh
        yh[ 2] = yoLEh
        yh[ 3] = ytLEh
        yh[ 4] = ytTEh
        yh[ 5] = yoTEh
        yh[ 6] = ycTEh
  
    # Vtail
    xv = zeros(6)
    yv = zeros(6)
    
    bov = parg[igbov]
    Sv  = parg[igSv]
    ARv = parg[igARv]
    lambdav = parg[iglambdav]
    sweepv  = parg[igsweepv]

    bv = parg[igbv]
    cov = parg[igcov]


    dx = parg[igxvbox]
    tanLv = tan(sweepv*π/180.0)
    ctv = cov*lambdav

    xaxv = 0.40

    xoLEv = cov*(    - xaxv) + dx
    xoTEv = cov*(1.0 - xaxv) + dx
    xtLEv = ctv*(    - xaxv) + dx + 0.5*(bv - bov)*tanLv
    xtTEv = ctv*(1.0 - xaxv) + dx + 0.5*(bv - bov)*tanLv

    yoLEv = 0.5*bov
    yoTEv = 0.5*bov
    ytLEv = 0.5*bv
    ytTEv = 0.5*bv

    xcLEv = xoLEv
    xcTEv = xoTEv
    ycLEv = yoLEv
    ycTEv = yoTEv

    xv[ 1] = xcLEv
    xv[ 2] = xoLEv
    xv[ 3] = xtLEv
    xv[ 4] = xtTEv
    xv[ 5] = xoTEv
    xv[ 6] = xcTEv

    yv[ 1] = ycLEv
    yv[ 2] = yoLEv
    yv[ 3] = ytLEv
    yv[ 4] = ytTEv
    yv[ 5] = yoTEv
    yv[ 6] = ycTEv

    # Fuel tank
        Rtank = Rfuse - 0.1 # Account for clearance_fuse
        l = parg[iglftankin]
        ARtank = 2.0
        xcyl0 = parg[igxftank] - l/2 + Rtank/ARtank
        xcyl1 = parg[igxftank] + l/2 - Rtank/ARtank
        ntank = 8
        xt = zeros(ntank*2 )
        yt = zeros(ntank*2 )
        for i = 1: ntank
            fraci = float(i-1)/float(ntank-1)
            fracx = cos(0.5*pi*fraci)

            k = i
            xt[k] = xcyl0 - Rtank/ARtank*fracx
            yt[k] = sqrt(Rtank^2 * max((1 - ((xt[k]-xcyl0)/(Rtank/ARtank))^2), 0.0) )
        end
        # k = k+1
        # xt[k] = xcyl0 + parg[iglftank]
        # yt[k] = Rtank
        for i = 1: ntank
            fraci = float(i-1)/float(ntank-1)
            fracx = sin(0.5*pi*fraci)

            k = i + ntank
            xt[k] = xcyl1 + (xcyl1 + Rtank/ARtank - xcyl1)*fracx
            yt[k] = sqrt(Rtank^2 * max((1 - ((xt[k]-xcyl1)/(Rtank/ARtank))^2), 0.0) )
        end

        # xt = LinRange(xcyl0 - Rfuse/ARtank , xcyl0, 20 )
        # yt = zeros(length(xt))
        # @. yt = sqrt(Rfuse^2 * max((1 - ((xt-xcyl0)/(Rfuse/ARtank))^2), 0.0) )

        xshell = zeros(ntank)
        yshell = zeros(ntank)
        AR = 3.0
        xshellcenter = parg[igxshell2] - Rfuse/AR
        for i = 1: ntank
            fraci = float(i-1)/float(ntank-1)
            fracx = sin(0.5*pi*fraci)

            k = i
            xshell[k] = xshellcenter + Rfuse/AR *fracx
            yshell[k] = sqrt(Rfuse^2 * max((1 - ((xshell[k]-xshellcenter)/(Rfuse/AR))^2), 0.0) )
        end

    #Seats
    pax = parg[igWpay]/parm[imWperpax]
    seat_pitch = 30.0 * in_to_m
    seat_width = 19.0 * in_to_m
    aisle_halfwidth = 10.0 * in_to_m # per CFR § 25.815 

    seats_per_row = Int(2*parg[igRfuse] ÷ (seat_width + aisle_halfwidth/3))
    rows = Int(ceil(pax / seats_per_row))
    
    println("Seats per row = $seats_per_row, Total rows = $rows")
    yseats, symmetric_seats = arrange_seats(seats_per_row, parg[igRfuse])
    
    if seats_per_row <= 10
        emergency_rows = [12, 13]
    else
        emergency_rows = [19, 20]
    end
    
    xseats = zeros(rows)'
    xseats[1] = parg[igxshell1 ] + 10.0*ft_to_m 
    for r in 2:rows
        emergency_exit = 0.0
        if (r in emergency_rows)
            emergency_exit = seat_pitch/2
        end
        xseats[r] = xseats[r-1] + seat_pitch + emergency_exit
    end

    ## Plot
    if ax === nothing
        # plt.style.use(["../miscellaneous/prash.mplstyle"]) # HACK
        fig, ax = plt.subplots(figsize=(8,5), dpi = 300)
    else
        ax.cla()
    end
    engz  = 25
    tankz = 10
    fusez = 5
    wingz = 6
    tailz = 1
        # Plot wing
            ax.plot(xw, yw, "-k", zorder = wingz)
            ax.plot(xw, -yw, "-k", zorder = wingz)
            
            # Panel break
            # ax.plot(xw[[2,5]],  yw[[2,5]], "-k", lw = 1.0, alpha = 0.5)
            # ax.plot(xw[[2,5]], -yw[[2,5]], "-k", lw = 1.0, alpha = 0.5)
        
        # Plot Tail
        parg[igzhtail] > 0 ? tailz = 21 : tailz = 1
            # ax.plot(xh,  yh, "-k", zorder = tailz)
            # ax.plot(xh, -yh, "-k", zorder = tailz)
            ax.fill_between(xh, -yh, yh, facecolor = "w", alpha = 0.8, edgecolor = "k", zorder = tailz, linewidth = 2.0)
        xvt = [-0.4, -0.3, -0.2, -0.15, 0.2, 0.6].*parg[igcov] .+ parg[igxvbox]
        tailthick = (parg[igcov]*parg[ighboxv]/2)
        yvt = hcat([0.0 0.5*tailthick 0.9*tailthick ones(2)' .*tailthick 0.0])[:]
        ax.fill_between(xvt, -yvt, yvt, facecolor = "k", alpha = 0.8, edgecolor = "k", zorder = 22)

        ax.plot(xv,yv, "--r", zorder = 21)

        # Plot fuse
            # ax.fill(xf,  yf, facecolor = "w", edgecolor = "k")
            # ax.fill(xf, -yf, facecolor = "w", edgecolor = "k")
            ax.fill_between(xf, -yf, yf, facecolor = "w", edgecolor = "k", zorder = fusez, linewidth = 2.0)
            
        # Tank
        if (pari[iifwing] == 0)
            ax.plot(xt,  yt, "k", lw = 1.5, alpha = 0.8, zorder = tankz)
            ax.plot(xt, -yt, "k", lw = 1.5, alpha = 0.8, zorder = tankz)
            ax.fill_between(xt, -yt, yt, facecolor = "r", alpha = 0.1, edgecolor = "k", zorder = tankz-1, linewidth = 1.0)
            ax.text(parg[igxftank], 0.0, "LH\$_2\$", fontsize = label_fs-2.0, zorder = tankz+1, ha="center", va="center")
        end

        # Xshell2
        ax.plot(xshell,  yshell, "k", lw = 1.5, zorder = tankz)
        ax.plot(xshell, -yshell, "k", lw = 1.5, zorder = tankz)

        # Plot Engines:
            D = parg[igdaftfan]
            
            lnace = parg[iglnaceaft]
            x = parg[igxtshaft] - 0.5
            ax.fill([x,x, x+lnace, x+lnace, x], [ D/8,  D/8 + D,  D/8 + D*3/4,  D/8 + 1/4*D,  D/8],
            lw = 1.5, edgecolor = "k", zorder = engz, facecolor = "w")
            ax.fill([x,x, x+lnace, x+lnace, x], [-D/8, -D/8 - D, -D/8 - D*3/4, -D/8 - 1/4*D, -D/8],
            lw = 1.5, edgecolor = "k", zorder = engz, facecolor = "w")

            ηs = bs/b
            ηo = bo/b
            D = parg[igdfan]
            neng = parg[igneng]
            lnace = parg[iglnace]
            dy = 2*D # space to leave near wing root and tip [m]
            if parg[igneng] == 2
                yi = [ηs*b/2]
            else
                yi = LinRange(bo/2 + dy , b/2 *3/4, Int(parg[igneng]/2))
            end
            xi = zero(yi)
            ηi = yi/(b/2)
            ci = zero(yi)
            for (i, η)  in enumerate(ηi)
                if η <=ηs
                    ci[i] = co*(1  + (λs -  1)*(η - ηo)/(ηs - ηo))
                else
                    ci[i] = co*(λs + (λt - λs)*(η - ηs)/(1  - ηs))
                end
            end

            tanL = tan(parg[igsweep]*π/180.0)
            @. xi = tanL * (yi - bo/2) - 0.4ci + parg[igxwbox] - 1.0

            ax.plot( [xi, xi, xi.+lnace, xi.+lnace, xi] ,      [yi.-D/2, yi.+D/2, yi.+D/3, yi.-D/3, yi.-D/2 ],
             color = "k", lw = 1.5, zorder = wingz-1)
            ax.plot( [xi, xi, xi.+lnace, xi.+lnace, xi] , -1 .*[yi.-D/2, yi.+D/2, yi.+D/3, yi.-D/3, yi.-D/2 ],
             color = "k", lw = 1.5, zorder = wingz-1)
            #Pylons
            
            ax.plot( [xi.+lnace/2, xi.+1.0] ,         [yi, yi], color = "k", lw = 2, zorder = wingz-2)
            ax.plot( [xi.+lnace/2, xi.+1.0] , -1.0 .* [yi, yi], color = "k", lw = 2, zorder = wingz-2)
      

        # Plot NP and CG range
            ax.scatter(parg[igxNP], 0.0, color = "k", marker="o", zorder = 21, label = "NP")
            ax.text(parg[igxNP], -1.0, "NP", fontsize=label_fs-2.0, ha="center", va="center", zorder = 21)

            ax.annotate("", xy=(parg[igxCGfwd ] , 0.0), xytext=(parg[igxCGaft ] , 0.0),
            fontsize=16, ha="center", va="bottom",
            arrowprops=Dict("arrowstyle"=> "|-|, widthA=0.2, widthB=0.2"),
            zorder = 21, label = "CG movement")
            ax.text(0.5*(parg[igxCGfwd ]+parg[igxCGaft ]), -1.0, "CG", fontsize=label_fs-2.0, ha="center", va="center", zorder = 21)

        # Show seats
        if symmetric_seats
            ax.scatter(ones(length(yseats),1).*xseats, ones(1,rows).* yseats, color = "gray", alpha = 0.1, marker = "s", s=15, zorder = 21)
            ax.scatter(ones(length(yseats),1).*xseats, ones(1,rows).*-yseats, color = "gray", alpha = 0.1, marker = "s", s=15, zorder = 21)
        else
            ax.scatter(ones(length(yseats),1).*xseats, ones(1,rows).* yseats, color = "gray", alpha = 0.1, marker = "s", s=15, zorder = 21)
        end
     # diagnostic marks
    #  ax.scatter(parg[igxftank] - l/2, 0.0, color = "k", marker="o", zorder = 21)
    #  ax.scatter(parg[igxftank], 0.0, color = "b", marker="o", zorder = 21)
    #  ax.scatter(parg[igxblend2], 0.0, color = "k", marker="o", zorder = 21)
    #  ax.plot([parg[igxftank]-l/2, parg[igxftank]+l/2],[0.0, 0.0], zorder = 21)


    # Annotations
    ax.text(0, 16, @sprintf("PFEI = %5.3f J/Nm\nM\$_{cruise}\$ = %.2f\nWMTO = %.1f tonnes\nSpan = %5.1f m\nco    = %5.1f m\n\$ \\Lambda \$ = %.1f\$^\\circ\$\nRfuse = %5.1f m\nL/D = %3.2f",
     parm[imPFEI], para[iaMach, ipcruise1],parg[igWMTO]/9.81/1000, parg[igb], parg[igco], parg[igsweep], parg[igRfuse], para[iaCL, ipcruise1]/para[iaCD, ipcruise1]),
     fontsize = label_fs, ha="left", va="top")

    yloc = -20
    ax.annotate("", xy=(0.0, yloc), xytext=( xf[end], yloc),
            fontsize=16, ha="center", va="bottom",
            arrowprops=Dict("arrowstyle"=> "|-|, widthA=0.5, widthB=0.5"),
             zorder = 30)
    ax.text(xend/2, yloc, @sprintf("l = %5.1f m", xend), bbox=Dict("ec"=>"w", "fc"=>"w"), ha="center", va="center", fontsize = 14, zorder = 31)

    # Span annotations:
    codeD = false
    codeE = false
    xcodeD = -2.0
    xcodeE = -3.5
        if codeD
            # ICAO code D 
            bmaxD = 36
            ax.vlines(xcodeD, -bmaxD/2, bmaxD/2, lw = 5, alpha = 0.2, color = "y")
            ax.hlines( bmaxD/2, xcodeD, 40.0, lw = 5, alpha = 0.2, color = "y")
            ax.hlines(-bmaxD/2, xcodeD, 40.0, lw = 5, alpha = 0.2, color = "y")
            ax.text(20, bmaxD/2+1, "ICAO Code D/ FAA Group III", color = "y", alpha = 0.8, fontsize = 12, ha="center", va="center")
        end
        if codeE
            # ICAO code E
            bmaxE = 52
            ax.vlines(xcodeE, -bmaxE/2, bmaxE/2, lw = 5, alpha = 0.2, color = "b")
            ax.hlines( bmaxE/2, xcodeE, 40.0, lw = 5, alpha = 0.2, color = "b")
            ax.hlines(-bmaxE/2, xcodeE, 40.0, lw = 5, alpha = 0.2, color = "b")
            ax.text(20, bmaxE/2+1, "ICAO Code E/ FAA Group IV", color = "b", alpha = 0.5, fontsize = 12, ha="center", va="center")
        end

    if codeE
        ax.set_ylim(-27,27)
    elseif codeD
        ax.set_ylim(-23,23)
    else
        ax.set_ylim(-23, 23)
    end
    ax.set_aspect(1)
    ax.set_xlim(-2.5, 52)
    ax.set_ylabel("y[m]")
    ax.set_xlabel("x[m]")
    plt.tight_layout()
    # ax.legend()
    ax.grid()

    if save_name !== nothing
        if pari[iifuel] == 1
            figname = @sprintf("ZIA_BLI_%d_%d_%.3f_%.1f", seats_per_row, parg[igneng], parm[imPFEI],  para[iaCL, ipcruise1]/para[iaCD, ipcruise1])
        elseif pari[iifuel] == 2
            figname = @sprintf("ZIA_SAF_BLI_%d_%d_%.3f_%.1f", seats_per_row, parg[igneng], parm[imPFEI],  para[iaCL, ipcruise1]/para[iaCD, ipcruise1])
        end
        plt.savefig(save_name*".png", metadata = Dict("Title"=>figname))
    end

    # Scale bar

    return ax
end

"""
    arrange_seats(seats_per_row, Rfuse,
     seat_width = 19.0 * in_to_m, 
     aisle_halfwidth = 10.0 * in_to_m)

Helper function to arrange seats given a number of `seats_per_row`
and fuselage radius. Assumes default `seat_width = 19"` and `aisle_halfwidth = 10"`,
but can be supplied by the user.
"""
function arrange_seats(seats_per_row, Rfuse,
     seat_width = 19.0 * in_to_m, 
     aisle_halfwidth = 10.0 * in_to_m)

    #Seats
    # Conditions:
    # - No more than 2 seats between any seat and the aisle
    seats_per_row % 2 == 0 ? symmetric_seats = true : symmetric_seats = false

    if symmetric_seats # seating can be symmetric
        half_seats_per_row = seats_per_row ÷ 2
        yseats = zeros(half_seats_per_row)

        if half_seats_per_row <= 3 #Single aisle
            yseats[1] = aisle_halfwidth + seat_width/2 #Aisle in the center
            for col = 2:half_seats_per_row
                yseats[col] = yseats[col-1] + seat_width #offset every seat by width
            end 
        else # twin aisle 
            #If symmetric no more than 2 seats next to each other at the 
            # centerline (I'm not evil enough to create a x-6-x seating arrangement even if "technically" allowed)
            yseats[1] = seat_width/2.0
            yseats[2] = yseats[1] + seat_width
            #Aisle
            half_seats_remaining = half_seats_per_row - 2
            if half_seats_remaining > 4
                @warn "Potentially trying to design a 3 aisle aircraft?
                Seating arrangement not (yet) automatically handled, so check carefully."
            end
            yseats[3] = yseats[2] + aisle_halfwidth*2 + seat_width
            for col = 4:half_seats_per_row
                yseats[col] = yseats[col-1] + seat_width
            end 
        end

    else
        @info "Asymmetric seating only deals with 3 or 5 seats at the moment"
        seating_excess_space = 2*Rfuse - seats_per_row*seat_width - 2*aisle_halfwidth 
        yseats = zeros(seats_per_row)
        # Start from edge and give some space based on the excess space available.
        ind = 1
        yseats[ind] = -Rfuse + seating_excess_space/2 + seat_width/2 
        ind+=1
        if seats_per_row > 3
            yseats[ind] = yseats[ind-1] + seat_width
            ind+=1
        end
        yseats[ind] = yseats[ind-1] + aisle_halfwidth*2 + seat_width
        ind+=1
        for col = ind:seats_per_row
            yseats[col] = yseats[col-1] + seat_width
        end 
    end
    return yseats, symmetric_seats
end  # function arrange_seats

"""
    PayloadRange(ac_og, Rpts, Ppts, filename, OEW, itermax)

Function to plot a payload range diagram for an aircraft

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac_og::aircraft`: Aircraft structure for payload range diagram.
    - `Rpts::Int64`: Density of ranges to be plot (Optional).
    - `Ppts::Int64`: Density of payloads to be plot (Optional).
    - `filename::String`: filename string for the plot to be stored (Optional).
    - `OEW::Boolean`: Whether to have OEW+Payload on the y-axis or just Payload (Optional).
    - `itermax::Int64`: Max Iterations for woper loop (Optional).
"""

function PayloadRange(ac_og, Rpts = 20, Ppts = 20, filename = "PayloadRangeDiagram.png", OEW = false, itermax = 20.0)
    ac = deepcopy(ac_og)
    RangeArray = ac.parm[imRange] * LinRange(0.1,1.2,Rpts)
    maxPay = 0

    Wmax = ac.parg[igWMTO]
    Fuelmax = ac.parg[igWfmax]
    Wempty = ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay]

    RangesToPlot = []
    PayloadToPlot = []
    maxPay = ac.parm[imWpay ]

    for Range = RangeArray
        Payloads = (maxPay) * LinRange(1, 0.1, Ppts)
        ac.parm[imRange] = Range
        for mWpay = Payloads
            println("Checking for Range (nmi): ",Range/1852.0, " and Pax = ", mWpay/(215*4.44822))
            ac.parm[imWpay ] = mWpay
            try
                @views TASOPT.woper(ac.pari,ac.parg,ac.parm[:,1:1],ac.para[:,:,1:1],ac.pare[:,:,1:1], ac.para[:,:,1:1],ac.pare[:,:,1:1], itermax,0.0)
                # woper success: store maxPay, break loop
                WTO = Wempty + mWpay + ac.parm[imWfuel]
                mWfuel = ac.parm[imWfuel]

                if WTO > Wmax || mWfuel > Fuelmax || WTO < 0.0 || mWfuel < 0.0 
                    WTO = 0.0
                    mWfuel = 0.0
                    println("Max out error!")
                else
                    maxPay = mWpay
                    println("Converged - moving to next range...")
                    break
                end     
            catch
                println("Not Converged - moving to lower payload...")      
            end
        end
        append!(RangesToPlot, Range)
        if OEW
            append!(PayloadToPlot, maxPay+Wempty)
        else
            append!(PayloadToPlot, maxPay)
        end
    end
    println(RangesToPlot)
    println(PayloadToPlot)
    fig, ax = plt.subplots(figsize=(8,5), dpi = 300)
    ax.plot(RangesToPlot ./ (1000*1852.0), PayloadToPlot./ (9.8*1000), linestyle="-",  color="b", label="Payload ")
    ax.set_xlabel("Range (1000 nmi)")
    ax.set_ylabel("Weight (1000 kg)")
    ax.legend()
    ax.set_title("Payload Range Plot")
    ax.grid()

    fig.savefig(filename)
end