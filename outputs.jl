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

function stickfig(parg, pari; ax = nothing)

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
        tanLh = tan(sweeph*pi/180.0)
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

        if  (true)
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
        xcyl0 = parg[igxftank] - parg[iglftank]/2
        xcyl1 = xcyl0 + parg[iglftank]
        ARtank = 2.0
        ntank = 5
        xt = zeros(ntank*2 )
        yt = zeros(ntank*2 )
        Rtank = Rfuse - 0.1 # Account for clearance_fuse
        for i = 1: ntank
            fraci = float(i-1)/float(ntank-1)
            fracx = cos(0.5*pi*fraci)

            k = i
            xt[k] = xcyl0 + (xcyl0 - Rtank/ARtank - xcyl0)*fracx
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


    ## Plot
    if ax === nothing
        plt.style.use(["~/prash.mplstyle"])
        fig, ax = plt.subplots(figsize=(8,5), dpi = 100)
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
            ax.plot(xh, yh, "-k")
            ax.plot(xh, -yh, "-k")

        # Plot fuse
            # ax.fill(xf, yf, facecolor = "w", edgecolor = "k")
            # ax.fill(xf, yf, facecolor = "w", edgecolor = "k")
            ax.fill_between(xf, -yf, yf, facecolor = "w", edgecolor = "k", zorder = 5, linewidth = 2.0)
            
            # Tank
            # ax.plot(xt,  yt, "k", lw = 1.5, zorder = 10)
            # ax.plot(xt, -yt, "k", lw = 1.5, zorder = 10)
            ax.fill_between(xt, -yt, yt, facecolor = "r", alpha = 0.1, edgecolor = "k", zorder = 6, linewidth = 1.0)
            ax.text(parg[igxftank], 0.0, "LH\$_2\$", fontsize = 14, zorder = 10, ha="center", va="center")

        # Plot NP and CG range
            ax.scatter(parg[igxNP], 0.0, color = "k", marker="o", zorder = 21, label = "NP")
            ax.text(parg[igxNP], -1.0, "NP", fontsize=14, ha="center", va="center", zorder = 21)

            ax.annotate("", xy=(parg[igxCGfwd ] , 0.0), xytext=(parg[igxCGaft ] , 0.0),
            fontsize=16, ha="center", va="bottom",
            arrowprops=Dict("arrowstyle"=> "|-|, widthA=0.2, widthB=0.2"),
            zorder = 21, label = "CG movement")
            ax.text(0.5*(parg[igxCGfwd ]+parg[igxCGaft ]), -1.0, "CG", fontsize=14, ha="center", va="center", zorder = 21)

    ax.set_aspect(1)
    ax.set_ylabel("y[m]")
    ax.set_xlabel("x[m]")

    # Annotations
    ax.text(0, 15, @sprintf("Span = %5.1f m\nco    = %5.1f m", parg[igb], parg[igco]),
     fontsize = 16, ha="left", va="top")

    yloc = -12
    ax.annotate("", xy=(0.0, yloc), xytext=( xf[end], yloc),
            fontsize=16, ha="center", va="bottom",
            arrowprops=Dict("arrowstyle"=> "|-|, widthA=0.5, widthB=0.5"),
             zorder = 30)
    ax.text(xend/2, yloc, @sprintf("l = %5.1f m", xend), bbox=Dict("ec"=>"w", "fc"=>"w"), ha="center", va="center", fontsize = 14, zorder = 31)
    plt.tight_layout()
    # ax.legend()

    return ax
end