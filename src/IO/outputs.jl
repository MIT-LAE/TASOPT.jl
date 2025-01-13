using PythonCall
"""
`weight_buildup` prints out the weight build up for the aircraft
"""
function weight_buildup(ac::aircraft; io=stdout)
    parg = ac.parg
    pari = ac.pari
    fuselage = ac.fuselage
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
    Whpesys = parg[igWMTO] * fuselage.HPE_sys.W
    Wlgnose = parg[igWMTO] * parg[igflgnose]
    Wlgmain = parg[igWMTO] * parg[igflgmain]
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

"""
    stickfig(parg,para,pari,parm; ax = nothing, label_fs = 16)

`stickfig` plots a "stick figure" airplane on the axis `ax` if provided 
(else a new axis is created). This is used in conjuction with other functions
like [`plot_details`](@ref) to create summary plots to track progress of optimization
or present results.
"""
function stickfig(ac::aircraft; ax = nothing, label_fs = 16, 
    annotate_text = true, annotate_length = true, annotate_group = true, show_grid = false)

    pari = ac.pari
    parg = ac.parg
    @views pare = ac.pare[:,:,1]
    @views para = ac.para[:,:,1]
    @views parm = ac.parm[:,:,1]
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    fuselage = ac.fuselage
    # Wing
        co = wing.layout.root_chord
        cs = wing.layout.root_chord*wing.inboard.λ
        ct = wing.layout.root_chord*wing.outboard.λ

        sweep = wing.layout.sweep
        λs = wing.inboard.λ
        λt = wing.outboard.λ

        bo = wing.layout.root_span
        bs = wing.layout.break_span
        b  = wing.layout.span

        xax = 0.40
        xcLE = -xax
        xcTE = 1.0 - xax

        dx = wing.layout.box_x
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
        
        Rfuse = fuselage.layout.radius
        wfb   = fuselage.layout.bubble_center_y_offset

        anose    = fuselage.layout.nose_radius
        btail    = fuselage.layout.tail_radius

        xnose    = fuselage.layout.x_nose
        xend     = fuselage.layout.x_end
        xblend1  = fuselage.layout.x_start_cylinder
        xblend2  = fuselage.layout.x_end_cylinder
        xhbox    = htail.layout.box_x
        
        hwidth = Rfuse + wfb
        
        nnose = 10
        ntail = 10

        xf = zeros(nnose + ntail + 1)
        yf = zeros(nnose + ntail + 1)

        if fuselage.layout.taper_fuse == 0
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
        
        boh = htail.layout.root_span
        Sh  = htail.layout.S
        ARh = htail.layout.AR
        lambdah = htail.outboard.λ
        sweeph  = htail.layout.sweep

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


        if (fuselage.layout.taper_fuse == 0)
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
  
        #Initialize seat start x-position
        xseats0 = fuselage.layout.x_pressure_shell_fwd
        # Fuel tank
        ntank = 8
        Rtank = Rfuse - 0.1 # Account for clearance_fuse
        l = max(parg[iglftankin], parg[iglftank])
        nftanks = pari[iinftanks] #Number of fuel tanks
        ARtank = 2.0

        if nftanks != 0
            tank_placement = ac.fuse_tank.placement
            if tank_placement == "front"
                xtanks = [parg[igxftank]]
                xseats0 = xtanks[1] + l/2 + 1.0 * ft_to_m #move seats backwards
            elseif tank_placement == "rear"
                xtanks = [parg[igxftankaft]]
                xseats0 = fuselage.layout.x_pressure_shell_fwd
            elseif tank_placement == "both"
                xtanks = [parg[igxftank], parg[igxftankaft]]
                xseats0 = xtanks[1] + l/2 + 1.0 * ft_to_m #move seats backwards
            end
            
            xt = zeros(nftanks, ntank*2 )
            yt = zeros(nftanks, ntank*2 )
            for m = 1:nftanks
                xcyl0 = xtanks[m] - l/2 + Rtank/ARtank
                xcyl1 = xtanks[m] + l/2 - Rtank/ARtank
                
                for i = 1: ntank
                    fraci = float(i-1)/float(ntank-1)
                    fracx = cos(0.5*pi*fraci)

                    k = i
                    xt[m, k] = xcyl0 - Rtank/ARtank*fracx
                    yt[m, k] = sqrt(Rtank^2 * max((1 - ((xt[m, k]-xcyl0)/(Rtank/ARtank))^2), 0.0) )
                end
                # k = k+1
                # xt[k] = xcyl0 + parg[iglftank]
                # yt[k] = Rtank
                for i = 1: ntank
                    fraci = float(i-1)/float(ntank-1)
                    fracx = sin(0.5*pi*fraci)

                    k = i + ntank
                    xt[m, k] = xcyl1 + (xcyl1 + Rtank/ARtank - xcyl1)*fracx
                    yt[m, k] = sqrt(Rtank^2 * max((1 - ((xt[m, k]-xcyl1)/(Rtank/ARtank))^2), 0.0) )
                end
            end
        end

        # xt = LinRange(xcyl0 - Rfuse/ARtank , xcyl0, 20 )
        # yt = zeros(length(xt))
        # @. yt = sqrt(Rfuse^2 * max((1 - ((xt-xcyl0)/(Rfuse/ARtank))^2), 0.0) )

        xshell = zeros(ntank)
        yshell = zeros(ntank)
        AR = 3.0
        xshellcenter = fuselage.layout.x_pressure_shell_aft - Rfuse/AR
        for i = 1: ntank
            fraci = float(i-1)/float(ntank-1)
            fracx = sin(0.5*pi*fraci)

            k = i
            xshell[k] = xshellcenter + Rfuse/AR *fracx
            yshell[k] = sqrt(Rfuse^2 * max((1 - ((xshell[k]-xshellcenter)/(Rfuse/AR))^2), 0.0) )
        end

    if pari[iidoubledeck] == 0 #Only show seats in single deck arrangements
        h_seat = fuselage.cabin.seat_height
        pax = parg[igWpay]/parm[imWperpax]
        Rfuse = fuselage.layout.radius
        dRfuse = fuselage.layout.bubble_lower_downward_shift
        wfb = fuselage.layout.bubble_center_y_offset
        nfweb = fuselage.layout.n_webs

        θ = find_floor_angles(false, Rfuse, dRfuse, h_seat = h_seat) #Find the floor angle
        wcabin = find_cabin_width(Rfuse, wfb, nfweb, θ, h_seat) #Cabin width
        _, xseats, seats_per_row = place_cabin_seats(pax, wcabin)

        xseats = xseats .+ xseats0
        rows = length(xseats)

        println("Seats per row = $seats_per_row, Total rows = $rows")
        yseats = arrange_seats(seats_per_row, wcabin)
    end
    ## Plot
    if ax === nothing
        # plt.style.use(["../miscellaneous/prash.mplstyle"]) # HACK
        fig, ax = plt.subplots(figsize=(8,5), dpi = 300)
        fig.subplots_adjust(right=0.75)
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
        htail.layout.z> 0 ? tailz = 21 : tailz = 1
            # ax.plot(xh,  yh, "-k", zorder = tailz)
            # ax.plot(xh, -yh, "-k", zorder = tailz)
            ax.fill_between(xh, -yh, yh, facecolor = "w", alpha = 0.8, edgecolor = "k", zorder = tailz, linewidth = 2.0)
        xvt = [-0.4, -0.15, 0.2, 0.6].*vtail.layout.root_chord .+ vtail.layout.box_x
        yvt = hcat([0.0 ones(length(xvt) - 2)' .*(vtail.layout.root_chord*vtail.outboard.cross_section.thickness_to_chord/2) 0.0])[:]
        ax.fill_between(xvt, -yvt, yvt, facecolor = "k", alpha = 0.8, edgecolor = "k", zorder = 22)

        # Plot fuse
            # ax.fill(xf,  yf, facecolor = "w", edgecolor = "k")
            # ax.fill(xf, -yf, facecolor = "w", edgecolor = "k")
            ax.fill_between(xf, -yf, yf, facecolor = "w", edgecolor = "k", zorder = 5, linewidth = 2.0)
            
        # Tank
        if (pari[iifwing] == 0)
            for m = 1:nftanks
                ax.plot(xt[m,:],  yt[m,:], "k", lw = 1.5, zorder = 10)
                ax.plot(xt[m,:], -yt[m,:], "k", lw = 1.5, zorder = 10)
                ax.fill_between(xt[m,:], -yt[m,:], yt[m,:], facecolor = "r", alpha = 0.1, edgecolor = "k", zorder = 6, linewidth = 1.0)

                if pari[iifuel] == 11
                    fuelname = "CH\$_4\$"
                elseif pari[iifuel] == 40
                    fuelname = "LH\$_2\$"
                end
                ax.text(xtanks[m], 0.0, fuelname, fontsize = label_fs-2.0, zorder = 10, ha="center", va="center")
            end
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
            ηs = wing.layout.ηs
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

            tanL = tan(wing.layout.sweep*π/180.0)
            @. xi = tanL * (yi - bo/2) - 0.4ci + wing.layout.box_x - 1.0
            xlocations = vec(hcat([xi, xi, xi.+lnace, xi.+lnace, xi]...)) #hcat is to avoid this being an array of arrays
            ylocations = vec(hcat([yi.-D/2, yi.+D/2, yi.+D/3, yi.-D/3, yi.-D/2]...))
            ax.plot(xlocations, ylocations, color = "r", linewidth = 1.5)

        # Plot NP and CG range
            ax.scatter(parg[igxNP], 0.0, color = "k", marker="o", zorder = 21, label = "NP")
            ax.text(parg[igxNP], -1.0, "NP", fontsize=label_fs-2.0, ha="center", va="center", zorder = 21)

            ax.annotate("", xy=(parg[igxCGfwd ] , 0.0), xytext=(parg[igxCGaft ] , 0.0),
            fontsize=16, ha="center", va="bottom",
            arrowprops=Dict("arrowstyle"=> "|-|, widthA=0.2, widthB=0.2"),
            zorder = 21, label = "CG movement")
            ax.text(0.5*(parg[igxCGfwd ]+parg[igxCGaft ]), -1.0, "CG", fontsize=label_fs-2.0, ha="center", va="center", zorder = 21)

        # Show seats

        if pari[iidoubledeck] == 0 #Show seats in single deck case
            ax.scatter(ones(length(yseats),1).*xseats, ones(1,rows).* yseats, color = "gray", alpha = 0.1, marker = "s", s=15, zorder = 21)
        end
     # diagnostic marks
    #  ax.scatter(parg[igxftank] - l/2, 0.0, color = "k", marker="o", zorder = 21)
    #  ax.scatter(parg[igxftank], 0.0, color = "b", marker="o", zorder = 21)
    #  ax.scatter(parg[igxblend2], 0.0, color = "k", marker="o", zorder = 21)
    #  ax.plot([parg[igxftank]-l/2, parg[igxftank]+l/2],[0.0, 0.0], zorder = 21)



    # Annotations
    if annotate_text
        ax.text(1.05, 0.75, transform=ax.transAxes, @sprintf("PFEI = %5.3f\nM\$_{cruise}\$ = %.2f\nWMTO = %.1f t\nSpan = %5.1f m\nco    = %5.1f m\n\$ \\Lambda \$ = %.1f\$^\\circ\$\nRfuse = %5.1f m\nL/D = %3.2f",
        parm[imPFEI], para[iaMach, ipcruise1],parg[igWMTO]/9.81/1000, wing.layout.span, wing.layout.root_chord, wing.layout.sweep, fuselage.layout.radius, para[iaCL, ipcruise1]/para[iaCD, ipcruise1]),
        fontsize = label_fs, ha="left", va="top")
    end
    if annotate_length
        yloc = -20
        ax.annotate("", xy=(0.0, yloc), xytext=( xf[end], yloc),
                fontsize=label_fs, ha="center", va="bottom",
                arrowprops=Dict("arrowstyle"=> "|-|, widthA=0.5, widthB=0.5"),
                zorder = 30)
        ax.text(xend/2, yloc, @sprintf("l = %5.1f m", xend), bbox=Dict("ec"=>"w", "fc"=>"w"), ha="center", va="center", fontsize = label_fs-2, zorder = 31)
    end
    # Span annotations:
    groups, bmax = find_aerodrome_code(wing.layout.max_span) #Find ICAO and FAA groups as well as max span
    xcode = -2.0

    if annotate_group
        ax.vlines(xcode, -bmax/2, bmax/2, lw = 5, alpha = 0.2, color = "y")
        ax.hlines( bmax/2, xcode, 40.0, lw = 5, alpha = 0.2, color = "y")
        ax.hlines(-bmax/2, xcode, 40.0, lw = 5, alpha = 0.2, color = "y")
        ax.text(20, bmax/2+1, "ICAO Code $(groups[1])/ FAA Group $(groups[2])", color = "y", alpha = 0.8, fontsize = 12, ha="center", va="center")
    end
    ax.set_ylim(-39, 39)
    ax.set_xlim(0, 80)
    # println("LIMS ARE ",-1.2*bmax/2, " ", 1.2*bmax/2)
    ax.set_aspect(1)
    ax.set_ylabel("y[m]")
    ax.set_xlabel("x[m]")
    plt.tight_layout()
    # ax.legend()

    if show_grid
        ax.grid()
    end

    return ax
end

"""
    find_aerodrome_code(b)

`find_aerodrome_code` finds the airport code corresponding to a given aircraft wingspan. It returns 
the codes in the ICAO and FAA formats, as well as the maximum wingspan for these codes.
"""
function find_aerodrome_code(b::Float64)
    max_bs = sort(collect(keys(aerodrome_codes)))
    idx = 1
    for cand_maxb in max_bs
        if b >= cand_maxb 
            idx += 1
        end
    end
    maxb = max_bs[idx]
    groups = aerodrome_codes[maxb]
    return groups, Float64(maxb)
end

"""
    plot_details(parg, pari, para, parm; ax = nothing)

`plot_details` combines a [`stickfig`](@ref) plot along with a mission summary,
weight and drag buildup stacked bar charts to present results.
"""
function plot_details(ac::aircraft; ax = nothing)

    pari = ac.pari
    parg = ac.parg
    fuselage = ac.fuselage
    wing = ac.wing
    htail = ac.wing
    vtail = ac.wing
    @views pare = ac.pare[:,:,1]
    @views para = ac.para[:,:,1]
    @views parm = ac.parm[:,:,1]
        ## Create empty plot
        if ax === nothing
            axd = plt.figure(figsize=(8,5), dpi = 300, layout="constrained").subplot_mosaic(
                    """
                    AA
                    BC
                    """,
                    # set the height ratios between the rows
                    height_ratios=[1, 3],
                    # set the width ratios between the columns
                    width_ratios=[1, 3],
                    )
            ax = [axd["A"], axd["B"], axd["C"]]
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
        Whpesys = parg[igWMTO] * fuselage.HPE_sys.W
        Wlgnose = parg[igWMTO] * parg[igflgnose]
        Wlgmain = parg[igWMTO] * parg[igflgmain]
        Wtotadd = Whpesys + Wlgnose + Wlgmain
        
        Wpay  = parg[igWpay]
        Wfuel = parg[igWfuel]
        WMTO  = parg[igWMTO]

        Wwing  = wing.weight
        Wfuse  = fuselage.weight
        Wvtail = vtail.weight
        Whtail = htail.weight
        Weng = parg[igWeng]
        Wtesys = parg[igWtesys]
        Wftank = parg[igWftank]

        Wwingfrac = Wwing /WMTO
        Wfusefrac = Wfuse /WMTO
        Wvtailfrac = Wvtail/WMTO
        Whtailfrac = Whtail/WMTO
        Wtesysfrac = Wtesys/WMTO
        Wftankfrac = Wftank/WMTO
        Wtotaddfrac = Wtotadd/WMTO
        Wengfrac = Weng/WMTO

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
        push!(Webars, a.bar(3, Wengfrac , bottom = Wfusefrac+Wwingfrac+Whtailfrac+Wvtailfrac+Wtesysfrac+Wftankfrac+Wtotaddfrac, width = bar_width, label = "Weng"))
        push!(Webars, a.bar(3, Wtotaddfrac , bottom = Wfusefrac+Wwingfrac+Whtailfrac+Wvtailfrac+Wtesysfrac+Wftankfrac, width = bar_width, label = "Wadd"))
        push!(Webars, a.bar(3, Wftankfrac  , bottom = Wfusefrac+Wwingfrac+Whtailfrac+Wvtailfrac+Wtesysfrac, width = bar_width, label = "Wftank"))
        push!(Webars, a.bar(3, Wtesysfrac  , bottom = Wfusefrac+Wwingfrac+Whtailfrac+Wvtailfrac, width = bar_width, label = "Wtesys"))
        push!(Webars, a.bar(3, Wvtailfrac  , bottom = Wfusefrac+Wwingfrac+Whtailfrac, width = bar_width, label = "Wvtail"))
        push!(Webars, a.bar(3, Whtailfrac  , bottom = Wfusefrac+Wwingfrac, width = bar_width, label = "Whtail"))
        push!(Webars, a.bar(3, Wwingfrac   , bottom = Wfusefrac, width = bar_width, label = "Wwing"))
        push!(Webars, a.bar(3, Wfusefrac   , width = bar_width, label = "Wfuse"))
        
        Welabels = ["Weng" "Wadd" "Wftank" "Wtesys" "Wvtail" "Whtail" "Wwing" "Wfuse"]
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
        w, h = bar[0].get_width(), bar[0].get_height()
        x, y = bar[0].get_x(), bar[0].get_y()
        a.text(x+w, y+h/2, @sprintf("%7.3f", pyconvert(Float64, h)*val_multiplier), ha = "left", va = "center", fontsize = fontsize)
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
    PayloadRange(ac_og; Rpts = 20, Ppts = 20, filename = "PayloadRangeDiagram.png", OEW = false, uniform_map = true)

Function to plot a payload range diagram for an aircraft

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac_og::aircraft`: Aircraft structure for payload range diagram.
    - `Rpts::Int64`: Density of ranges to be plot (Optional).
    - `Ppts::Int64`: Density of payloads to be plot (Optional).
    - `filename::String`: filename string for the plot to be stored (Optional).
    - `OEW::Boolean`: Whether to have OEW+Payload on the y-axis or just Payload (Optional).
    - `uniform_map::Boolean`: use a uniform or non-uniform spacing between ranges and payloads (Optional)
"""
function PayloadRange(ac_og; Rpts = 20, Ppts = 20, filename = "PayloadRangeDiagram.png", OEW = false, uniform_map = true)
    #Duplicate design mission as second mission, which will be modified
    parm = cat(ac_og.parm[:,1], ac_og.parm[:,1], dims=2)
    pare = cat(ac_og.pare[:,:,1], ac_og.pare[:,:,1], dims=3)
    para = cat(ac_og.para[:,:,1], ac_og.para[:,:,1], dims=3)
    ac = aircraft(ac_og.name, ac_og.description,
    ac_og.pari, ac_og.parg, parm, para, pare, [true], ac_og.fuse_tank, ac_og.fuselage)

    #Extract aircraft parameters
    maxPay = ac.parg[igWpaymax]
    Wmax = ac.parg[igWMTO]
    Fuelmax = ac.parg[igWfmax]
    Wempty = ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay]
    PFEI = 0

    #Use a uniform or non-uniform map for the range 
    if uniform_map #Use unifrom map
        RangeMap = LinRange(0.1,2, Rpts)
        PayMap = LinRange(1, 0, Ppts)
    else #greater density around design range
        RangeMap = unique([LinRange(0.01,1, Int64(ceil(Rpts/2)+1)).^0.5; 1 .+ LinRange(0,1, Int64(floor(Rpts/2))).^2])
        PayMap = LinRange(1, 0, Ppts).^0.75 #Use a map with greater density in the high payloads
    end
    RangeArray = ac.parm[imRange,1] * RangeMap
    Payloads = PayMap * maxPay

    RangesToPlot = []
    PayloadToPlot = []
    PFEIsToPlot = []

    for Range = RangeArray
        if maxPay == 0
            break   
        end
        ac.parm[imRange,2] = Range
        for mWpay = Payloads
            println("Checking for Range (nmi): ",Range/1852.0, " and Pax = ", mWpay/(215*4.44822))
            ac.parm[imWpay,2] = mWpay
            try
                TASOPT.fly_off_design(ac, 2)
                # fly_off_design success: store maxPay, break loop
                mWfuel = ac.parm[imWfuel,2]
                WTO = Wempty + mWpay + mWfuel

                if WTO > Wmax || mWfuel > Fuelmax || WTO < 0.0 || mWfuel < 0.0 
                    WTO = 0.0
                    mWfuel = 0.0
                    println("Max out error!")
                    if mWpay == 0
                        println("Payload 0 and no convergence found")
                        maxPay = 0
                        PFEI = NaN
                    end
                else
                    maxPay = mWpay
                    PFEI = ac.parm[imPFEI, 2]
                    println("Converged - moving to next range...")
                    break
                end     
            catch
                println("Not Converged - moving to lower payload...")      
            end
        end
        append!(RangesToPlot, Range)
        append!(PFEIsToPlot, PFEI)
        if OEW
            append!(PayloadToPlot, maxPay+Wempty)
        else
            append!(PayloadToPlot, maxPay)
        end
    end
    println(RangesToPlot)
    println(PayloadToPlot)
    fig, axes = pyplot.subplots(2,1,figsize=(8,10), dpi = 300)
    ax1 = axes[0]
    ax2 = axes[1]

    ax1.plot(RangesToPlot ./ (1000*1852.0), PayloadToPlot./ (9.8*1000), linestyle="-",  color="b", label="Payload ")
    ax1.set_xlabel("Range (1000 nmi)")
    ax1.set_ylabel("Weight (1000 kg)")
    ax1.legend()
    ax1.set_xlim([0,RangesToPlot[end] / (1000*1852.0)])
    ax1.set_title("Payload-Range Plot")
    ax1.grid()

    ax2.plot(RangesToPlot ./ (1000*1852.0), PFEIsToPlot, linestyle="-",  color="b")
    ax2.set_xlabel("Range (1000 nmi)")
    ax2.set_ylabel("PFEI at maximum payload")
    ax2.legend()
    ax2.set_ylim([0,2])
    ax2.set_xlim([0,RangesToPlot[end] / (1000*1852.0)])
    ax2.set_title("PFEI-Range Plot")
    ax2.grid()

    fig.savefig(filename)
end