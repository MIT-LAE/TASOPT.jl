"""
    size_wing_section!(layout, section, cs, cp, 
        cap_material, tauweb, sigfac)

Calculates Loads and thicknesses for wing sections

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `layout::TASOPT.structures.Wing.WingLayout`: Wing Layout.
    - `section::TASOPT.structures.Wing.WingSection`: Wing Section to be sized.
    - `cs::Float64`: Wing section chord.
    - `cp::Float64`: Chord times cosine of sweep
    - `cap_material::TASOPT.materials`: Material of cap.
    - `tauweb::Float64`: Webs tau
    - `sigfac::Float64`: Stress factor
"""
function size_wing_section!(layout, section, cs, cp, cap_material, tauweb, sigfac)
    shear_load = section.max_shear_load
    moment = section.moment

    Eweb = cap_material.E
    Ecap = cap_material.E
    Gcap = cap_material.G
    Gweb = cap_material.G
    sigcap = cap_material.σmax * sigfac

    cosL = cos(layout.sweep*pi/180)

    h_avg, h_rms = get_average_sparbox_heights(section.layout)

    tbweb = calc_web_thickness(tauweb, shear_load, cs*cosL, layout.hweb_to_hbox*section.layout.thickness_to_chord)
    tbcap = calc_cap_thickness(sigcap, moment, section.layout.thickness_to_chord,
        layout.box_width, h_rms, cs, cosL)

    Abcap = 2.0*tbcap*layout.box_width
    Abweb = 2.0*tbweb*layout.hweb_to_hbox*section.layout.thickness_to_chord

    section.web_cap.EI_bending = Ecap*cp^4 * (h_rms^3 - (h_rms-2.0*tbcap)^3)*layout.box_width/12.0 +
        Eweb*cp^4 * tbweb*(layout.hweb_to_hbox*section.layout.thickness_to_chord)^3 / 6.0
    
    section.web_cap.EI_normal = Ecap*cp^4 * tbcap*    layout.box_width  ^3 / 6.0 +
        Eweb*cp^4 * tbweb*layout.hweb_to_hbox*section.layout.thickness_to_chord * 0.5*layout.box_width^2
            
    section.web_cap.GJ = cp^4 * 2.0*((layout.box_width-tbweb)*(h_avg-tbcap))^2 /
        (  (layout.hweb_to_hbox*section.layout.thickness_to_chord-tbcap)/(Gweb*tbweb) +
        (   layout.box_width -tbweb)/(Gcap*tbcap) )

    return tbweb, tbcap, Abcap, Abweb  
end

"""
    surfw!(wing, po, gammat, gammas, 
    Nload, We, neout, dyeout, neinn, dyeinn,
    fLt, sigfac, rhofuel)

Calculates Wing or Tail loads, stresses, weights of individual wing sections.
Also returns the material gauges, torsional and bending stiffness.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure.
    - `po::Float64`: Point where loads and stresses are calculated.
    - `gammat::Float64`: Tip airfoil section shape exponent.
    - `gammas::Float64`: Start airfoil section shape exponent.
    - `Nload::Int`: Number of loads (used to distribute engine loads).
    - `iwplan::Int`: Indicates the presence of a strut.
    - `We::Float64`: Weight of the engine.
    - `neout::Int`: Number of outboard engines.
    - `dyeout::Float64`: Distance between engines and the wingtip.
    - `neinn::Int`: Number of inboard engines.
    - `dyeinn::Float64`: Distance between engines and the wing root.
    - `fLt::Float64`: Factor applied to the tip load.
    - `sigfac::Float64`: Stress Factor.
    - `rhofuel::Float64`: Density of the fuel.

See [Geometry](@ref geometry),  [Wing/Tail Structures](@ref wingtail), and Section 2.7  of the [TASOPT Technical Description](@ref dreladocs). 
"""
function surfw!(wing, po, gammat, gammas, 
       Nload, We, neout, dyeout, neinn, dyeinn,
       fLt, sigfac, rhofuel)

    tauweb,sigstrut = wing.inboard.webs.material.τmax * sigfac, wing.strut.material.σmax * sigfac

    cosL = cos(wing.layout.sweep*pi/180)
    sinL = sin(wing.layout.sweep*pi/180)

    # Calculate non-dim span coordinate at span break and root (ηs and ηo resp.)
    etao = wing.outboard.layout.b/wing.layout.b
    etas = wing.inboard.layout.b/wing.layout.b

    cop = wing.layout.chord*cosL
    csp = wing.layout.chord*cosL*wing.inboard.layout.λ
    # Tip roll off Lift (modeled as a point load) and it's moment about ηs
    dLt = fLt*po*wing.layout.chord*gammat*wing.outboard.layout.λ
    dMt = dLt*0.5*wing.layout.b*(1.0-etas)

    h_avgo = wing.inboard.layout.thickness_to_chord * (1.0 - (1.0-wing.layout.hweb_to_hbox)/3.0)
    h_avgs = wing.outboard.layout.thickness_to_chord * (1.0 - (1.0-wing.layout.hweb_to_hbox)/3.0)

    # Outboard section:
    #---- strut-attach shear,moment from outer-wing loading. 
    #     Note added term to account for any outboard engines.
    #     If neout = 0 this simplifies to Drela's version which assumes engine
    #     fixed at ηs locations only.
    wing.outboard.max_shear_load = (po*wing.layout.b   / 4.0)*(gammas+    gammat)*(1.0-etas) +
    dLt - Nload*wing.outboard.weight - Nload*neout*We
    wing.outboard.moment = (po*wing.layout.b^2/24.0)*(gammas+2.0*gammat)*(1.0-etas)^2 +
    dMt - Nload*wing.outboard.dyW - Nload*neout*We*dyeout

    #---- size strut-attach station at etas
    cs = wing.layout.chord*wing.inboard.layout.λ

    tbwebs, tbcaps, Abcaps, Abwebs = size_wing_section!(wing.layout, wing.outboard, cs,csp, wing.inboard.caps.material, tauweb, sigfac)
    # Inboard Section:
    if(wing.planform==0 || wing.planform==1) 
        #----- no strut, with or without engine at etas
        ls = 0.
        Tstrut = 0.
        Rstrut = 0.
        Pstrut = 0.

        # Modifed to account for bending relief from multiple engines.
        # dyeinn allows engine to be in locations other than ηs
        So = wing.outboard.max_shear_load - Nload*neinn*We +
            0.25*po*wing.layout.b*(1.0+gammas)*(etas-etao) -
            Nload*wing.inboard.weight
        Mo = wing.outboard.moment + wing.outboard.max_shear_load*0.5*wing.layout.b*(etas-etao) +
            (1.0/24.0)*po*wing.layout.b^2*(1.0+2.0*gammas)*(etas-etao)^2 -
            Nload*wing.inboard.dyW - Nload*neinn*We*dyeinn

        #----- limit So,Mo to Ss,Ms, which might be needed with heavy outboard engine
        #-      (rules out negatively-tapered structure, deemed not feasible for downloads)
        wing.inboard.max_shear_load = max(So ,wing.outboard.max_shear_load)
        wing.inboard.moment = max(Mo ,wing.outboard.moment)

        tbwebo, tbcapo, Abcapo, Abwebo = size_wing_section!(wing.layout, wing.inboard, wing.layout.chord, cop, wing.inboard.caps.material, tauweb, sigfac)

        lsp = 0.
        Tstrutp = 0.
        wing.strut.cos_lambda = 1.0

    else
        #----- strut present
        ls = sqrt(wing.strut.z^2 + (0.5*wing.layout.b*(etas-etao))^2)
        Rstrut = (po*wing.layout.b/12.0)*(etas-etao)*(1.0+2.0*gammas) + wing.outboard.max_shear_load
        Tstrut = Rstrut*ls/wing.strut.z
        #c     Pstrut = Rstrut*0.5*wing.layout.b*(etas-etao)/zs

        #----- inboard shear,moment used for sparbox sizing
        wing.inboard.max_shear_load = wing.outboard.max_shear_load
        wing.inboard.moment = wing.outboard.moment
        #
        #----- size inboard station at etao
        tbwebo, tbcapo, Abcapo, Abwebo = size_wing_section!(wing.layout, wing.inboard, wing.layout.chord, cop, wing.inboard.caps.material, tauweb, sigfac)

        #----- total strut length, tension
        lsp = sqrt(wing.strut.z^2 + (0.5*wing.layout.b*(etas-etao)/cosL)^2)
        Tstrutp = Tstrut*lsp/ls
        wing.strut.cos_lambda = ls/lsp

    end

    Abfuels = (wing.layout.box_width-2.0*tbwebs)*(h_avgs-2.0*tbcaps)
    Abfuelo = (wing.layout.box_width-2.0*tbwebo)*(h_avgo-2.0*tbcapo)

    wing.strut.axial_force = Tstrutp/sigstrut

    Vcen = wing.layout.chord^2*wing.layout.b*  etao / 2.0

    Vinn = wing.layout.chord^2*wing.layout.b* (etas-etao) *
        (1.0 + wing.inboard.layout.λ + wing.inboard.layout.λ^2)/6.0 *
        cosL
    Vout = wing.layout.chord^2*wing.layout.b* (1.0 -etas) *
        (wing.inboard.layout.λ^2 + wing.inboard.layout.λ*wing.outboard.layout.λ + wing.outboard.layout.λ^2)/6.0 *
        cosL

    dxVinn = wing.layout.chord^2*wing.layout.b^2 * (etas-etao)^2 *
        (1.0 + 2.0*wing.inboard.layout.λ + 3.0*wing.inboard.layout.λ^2)/48.0 *
        sinL
    dxVout = wing.layout.chord^2*wing.layout.b^2 * (1.0 -etas)^2 *
        (wing.inboard.layout.λ^2 + 2.0*wing.inboard.layout.λ*wing.outboard.layout.λ + 3.0*wing.outboard.layout.λ^2)/48.0 *
        sinL +
        wing.layout.chord^2*wing.layout.b^2 * (etas-etao)*(1.0 -etas) *
        (wing.inboard.layout.λ^2 + wing.inboard.layout.λ*wing.outboard.layout.λ + wing.outboard.layout.λ^2)/12.0 *
        sinL

    dyVinn = wing.layout.chord^2*wing.layout.b^2 * (etas-etao)^2 *
        (1.0 + 2.0*wing.inboard.layout.λ + 3.0*wing.inboard.layout.λ^2)/48.0 *
        cosL
    dyVout = wing.layout.chord^2*wing.layout.b^2 * (1.0 -etas)^2 *
        (wing.inboard.layout.λ^2 + 2.0*wing.inboard.layout.λ*wing.outboard.layout.λ + 3.0*wing.outboard.layout.λ^2)/48.0 *
        cosL

    #---- set chord^2 weighted average areas for inner panel
    Abcapi = (Abcapo + Abcaps*wing.inboard.layout.λ^2)/(1.0+wing.inboard.layout.λ^2)
    Abwebi = (Abwebo + Abwebs*wing.inboard.layout.λ^2)/(1.0+wing.inboard.layout.λ^2)

    Wscen   = (wing.inboard.caps.ρ*Abcapo + wing.inboard.webs.ρ*Abwebo)*gee*Vcen
    Wsinn   = (wing.inboard.caps.ρ*Abcapi + wing.inboard.webs.ρ*Abwebi)*gee*Vinn
    Wsout   = (wing.inboard.caps.ρ*Abcaps + wing.inboard.webs.ρ*Abwebs)*gee*Vout

    dxWsinn = (wing.inboard.caps.ρ*Abcapi + wing.inboard.webs.ρ*Abwebi)*gee*dxVinn
    dxWsout = (wing.inboard.caps.ρ*Abcaps + wing.inboard.webs.ρ*Abwebs)*gee*dxVout

    dyWsinn = (wing.inboard.caps.ρ*Abcapi + wing.inboard.webs.ρ*Abwebi)*gee*dyVinn
    dyWsout = (wing.inboard.caps.ρ*Abcaps + wing.inboard.webs.ρ*Abwebs)*gee*dyVout


    Abfueli= (Abfuelo + Abfuels*wing.inboard.layout.λ^2)/(1.0+wing.inboard.layout.λ^2)

    Wfcen   = rhofuel*Abfuelo *gee*Vcen
    Wfinn   = rhofuel*Abfueli *gee*Vinn
    Wfout   = rhofuel*Abfuels *gee*Vout

    dxWfinn = rhofuel*Abfueli *gee*dxVinn
    dxWfout = rhofuel*Abfuels *gee*dxVout

    dyWfinn = rhofuel*Abfueli *gee*dyVinn
    dyWfout = rhofuel*Abfuels *gee*dyVout

    wing.inboard.caps.weight.W   = 2.0*wing.inboard.caps.ρ*gee*(Abcapo*Vcen + Abcapi*Vinn + Abcaps*Vout)
    wing.inboard.webs.weight.W   = 2.0*wing.inboard.webs.ρ*gee*(Abwebo*Vcen + Abwebi*Vinn + Abwebs*Vout)

    dxWcap = 2.0*wing.inboard.caps.ρ*gee*( Abcapi*dxVinn + Abcaps*dxVout )
    dxWweb = 2.0*wing.inboard.webs.ρ*gee*( Abwebi*dxVinn + Abwebs*dxVout )

    wing.strut.weight   = 2.0*wing.strut.material.ρ*gee*wing.strut.axial_force*lsp
    wing.strut.dxW = wing.strut.weight * 0.25*wing.layout.b*(etas-etao) * sinL/cosL

    Vout = wing.layout.chord^2*wing.layout.b* (1.0 -etas) *
    (wing.inboard.layout.λ^2 + wing.inboard.layout.λ*wing.outboard.layout.λ + wing.outboard.layout.λ^2)/6.0 *
    cosL

    wing.inboard.caps.thickness = tbcapo
    wing.inboard.webs.thickness = tbwebo
    wing.outboard.caps.thickness = tbcaps
    wing.outboard.webs.thickness = tbwebs

    fwadd = wing_additional_weight(wing)
    Wwing = 2.0 * (Wscen + Wsinn + Wsout) * (1.0 + fwadd)
    wing.dxW = 2.0 * (dxWsinn + dxWsout) * (1.0 + fwadd)

    return Wwing,Wsinn,Wsout,dyWsinn,dyWsout,Wfcen,Wfinn,Wfout,dxWfinn,dxWfout,dyWfinn,dyWfout,lsp

end # surfw


"""
"""
function calc_web_thickness(τmax, shear, c_perp, web_height)
    t_web = shear / (c_perp^2 * 2 * web_height * τmax)
    return t_web
end  # function calc_web_thickness

"""
"""
function calc_cap_thickness(σmax, moment, h̄, w̄, h_rms, c, cosΛ)
    con = moment * 6h̄ / w̄ * 1 / (c^3 * σmax * cosΛ^4)
    t_cap = 0.5 * (h_rms - ∛(h_rms^3 - con))
    return t_cap
end  # function calc_cap_thickness