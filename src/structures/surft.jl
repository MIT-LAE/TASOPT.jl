"""
    surft!(tail, po, lambdas,gammat,gammas,fLt,
            tauweb,sigcap,Ecap,rhoweb,rhocap, b = tail.layout.b)

Calculates Tail loads, stresses, weights of individual tail sections.
Also returns the material gauges, torsional and bending stiffness.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `tail::TASOPT.structures.htail`: Tail structure.
    - `poh::Float64`: Point where loads and stresses are calculated.
    - `lambdas::Float64`: Start chord ratio (start chord / root chord).
    - `gammat::Float64`: Tip airfoil section shape exponent.
    - `gammas::Float64`: Start airfoil section shape exponent.
    - `fLt::Float64`: Factor applied to the tip load.
    - `tauweb::Float64`: Web material shear strength.
    - `sigcap::Float64`: Cap material axial compressive strength.
    - `Ecap::Float64`: Cap material Young's modulus.
    - `rhoweb::Float64`: Density of the web material.
    - `rhocap::Float64`: Density of the cap material.
    - `rhofuel::Float64`: Density of the fuel.
    - `b::Float64`: Wingspan.

See [Geometry](@ref geometry),  [Wing/Tail Structures](@ref wingtail), and Section 2.7  of the [TASOPT Technical Description](@ref dreladocs). 
"""

function surft!(tail, po, lambdas,gammat,gammas,fLt,
            tauweb,sigcap,Ecap,rhoweb,rhocap, b = tail.layout.b)

       Eweb = Ecap
       Gcap = Ecap * 0.5 / (1.0 + 0.3)
       Gweb = Ecap * 0.5 / (1.0 + 0.3)

       cosL = cos(deg2rad(tail.layout.sweep))
       sinL = sin(deg2rad(tail.layout.sweep))

      # Calculate non-dim span coordinate at span break and root (ηs and ηo resp.)
      etao = tail.outboard.layout.b/b
      etas = tail.outboard.layout.b/b

      cop = tail.layout.chord*cosL
      # Tip roll off Lift (modeled as a point load) and it's moment about ηs
      dLt = fLt*po*tail.layout.chord*gammat*tail.outboard.layout.λ
      dMt = dLt*0.5*b*(1.0-etas)

      havgo = tail.outboard.layout.chord_thickness * (1.0 - (1.0-tail.layout.hweb_to_hbox)/3.0)

      hrmso = tail.outboard.layout.chord_thickness * sqrt(1.0 - (1.0-tail.layout.hweb_to_hbox)/1.5 + (1.0-tail.layout.hweb_to_hbox)^2 / 5.0)
      hrmss = tail.outboard.layout.chord_thickness * sqrt(1.0 - (1.0-tail.layout.hweb_to_hbox)/1.5 + (1.0-tail.layout.hweb_to_hbox)^2 / 5.0)

# Outboard section:
#---- strut-attach shear,moment from outer-wing loading. 
#     Note added term to account for any outboard engines.
#     If neout = 0 this simplifies to Drela's version which assumes engine
#     fixed at ηs locations only.
      Ss = (po*b   / 4.0)*(gammas+    gammat)*(1.0-etas) +
	 dLt 
      Ms = (po*b^2/24.0)*(gammas+2.0*gammat)*(1.0-etas)^2 +
	 dMt

#---- size strut-attach station at etas
      cs = tail.layout.chord*lambdas
      tbwebs = Ss*0.5/(cs^2*tauweb*tail.layout.hweb_to_hbox*tail.outboard.layout.chord_thickness*cosL^2)
      con  = Ms*6.0*tail.outboard.layout.chord_thickness/(cs^3*sigcap*tail.layout.box_width*cosL^4)
      tbcaps = 0.5*(hrmss - (hrmss^3-con)^(1.0/3.0))
      Abcaps = 2.0*tbcaps*tail.layout.box_width
      Abwebs = 2.0*tbwebs*tail.layout.hweb_to_hbox*tail.outboard.layout.chord_thickness

    #----- no strut, with or without engine at etas

    # Modifed to account for bending relief from multiple engines.
    # dyeinn allows engine to be in locations other than ηs
    So = Ss +
    0.25*po*b*(1.0+gammas)*(etas-etao) 
    Mo = Ms + Ss*0.5*b*(etas-etao) +
    (1.0/24.0)*po*b^2*(1.0+2.0*gammas)*(etas-etao)^2 

    #----- limit So,Mo to Ss,Ms, which might be needed with heavy outboard engine
    #-      (rules out negatively-tapered structure, deemed not feasible for downloads)
    So = max(So ,Ss)
    Mo = max(Mo ,Ms)

    #----- size root station at etao
    tbwebo = So*0.5/(tail.layout.chord^2*tauweb*tail.layout.hweb_to_hbox*tail.outboard.layout.chord_thickness*cosL^2)
    con  = Mo*6.0*tail.outboard.layout.chord_thickness/(tail.layout.chord^3*sigcap*tail.layout.box_width*cosL^4)
    tbcapo = 0.5*(hrmso - (hrmso^3-con)^(1.0/3.0))
    Abcapo = 2.0*tbcapo*tail.layout.box_width
    Abwebo = 2.0*tbwebo*tail.layout.hweb_to_hbox*tail.outboard.layout.chord_thickness


      EIco = Ecap*cop^4 * (hrmso^3 - (hrmso-2.0*tbcapo)^3)*tail.layout.box_width/12.0 +
	 Eweb*cop^4 * tbwebo*(tail.layout.hweb_to_hbox*tail.outboard.layout.chord_thickness)^3 / 6.0
      EIno = Ecap*cop^4 * tbcapo*    tail.layout.box_width  ^3 / 6.0 +
	 Eweb*cop^4 * tbwebo*tail.layout.hweb_to_hbox*tail.outboard.layout.chord_thickness * 0.5*tail.layout.box_width^2
      GJo = cop^4 * 2.0*((tail.layout.box_width-tbwebo)*(havgo-tbcapo))^2 /
	 (  (tail.layout.hweb_to_hbox*tail.outboard.layout.chord_thickness-tbcapo)/(Gweb*tbwebo) +
	 (   tail.layout.box_width -tbwebo)/(Gcap*tbcapo) )


      Vcen = tail.layout.chord^2*b *  etao / 2.0

      Vinn = tail.layout.chord^2*b * (etas-etao) *
	 (1.0 + lambdas + lambdas^2)/6.0 *
	 cosL
      Vout = tail.layout.chord^2*b * (1.0 -etas) *
	 (lambdas^2 + lambdas*tail.outboard.layout.λ + tail.outboard.layout.λ^2)/6.0 *
	 cosL

      dxVinn = tail.layout.chord^2*b^2 * (etas-etao)^2 *
	 (1.0 + 2.0*lambdas + 3.0*lambdas^2)/48.0 *
	 sinL
      dxVout = tail.layout.chord^2*b^2 * (1.0 -etas)^2 *
	 (lambdas^2 + 2.0*lambdas*tail.outboard.layout.λ + 3.0*tail.outboard.layout.λ^2)/48.0 *
	 sinL +
	  tail.layout.chord^2*b^2 * (etas-etao)*(1.0 -etas) *
	 (lambdas^2 + lambdas*tail.outboard.layout.λ + tail.outboard.layout.λ^2)/12.0 *
	 sinL

#---- set chord^2 weighted average areas for inner panel
      Abcapi = (Abcapo + Abcaps*lambdas^2)/(1.0+lambdas^2)
      Abwebi = (Abwebo + Abwebs*lambdas^2)/(1.0+lambdas^2)

      Wscen   = (rhocap*Abcapo + rhoweb*Abwebo)*gee*Vcen
      Wsinn   = (rhocap*Abcapi + rhoweb*Abwebi)*gee*Vinn
      Wsout   = (rhocap*Abcaps + rhoweb*Abwebs)*gee*Vout

      dxWsinn = (rhocap*Abcapi + rhoweb*Abwebi)*gee*dxVinn
      dxWsout = (rhocap*Abcaps + rhoweb*Abwebs)*gee*dxVout

      Vout = tail.layout.chord^2*b * (1.0 -etas) *
	 (lambdas^2 + lambdas*tail.outboard.layout.λ + tail.outboard.layout.λ^2)/6.0 *
	 cosL

        Whtail = tail.ntails * (Wscen + Wsinn + Wsout) * (1.0 + tail.weight_fraction_added)
        dxWhtail = tail.ntails * (dxWsinn + dxWsout) * (1.0 + tail.weight_fraction_added)
        tail.weight = Whtail
        tail.outboard.dxW = dxWhtail

     tail.outboard.thickness_web, tail.outboard.thickness_cap, tail.outboard.EI_bending, tail.outboard.EI_normal, tail.outboard.GJ =  tbwebo,tbcapo,EIco,EIno,GJo
end # surft!



