"""
    surfw(po, b, bs, bo, co, zs,
        tail.outboard.λ, lambdas, gammat, gammas,
        Nload, iwplan, We, neout, dyeout, neinn, dyeinn,
        Winn, Wout, dyWinn, dyWout,
        sweep, wbox, tail.outboard.chord_thickness, tail.outboard.chord_thickness, rh, fLt,
        tauweb, sigcap, sigstrut, Ecap, Eweb, Gcap, Gweb,
        rhoweb, rhocap, rhostrut, rhofuel)

Calculates Wing or Tail loads, stresses, weights of individual wing sections.
Also returns the material gauges, torsional and bending stiffness.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `po::Float64`: Point where loads and stresses are calculated.
    - `b::Float64`: Wingspan.
    - `bs::Float64`: Spanwise location of the start of the taper.
    - `bo::Float64`: Spanwise location of the root chord.
    - `co::Float64`: Root chord length.
    - `zs::Float64`: Height of the strut attach point above wing.
    - `tail.outboard.λ::Float64`: Tip chord ratio (tip chord / root chord).
    - `lambdas::Float64`: Start chord ratio (start chord / root chord).
    - `gammat::Float64`: Tip airfoil section shape exponent.
    - `gammas::Float64`: Start airfoil section shape exponent.
    - `Nload::Int`: Number of loads (used to distribute engine loads).
    - `iwplan::Int`: Indicates the presence of a strut.
    - `We::Float64`: Weight of the engine.
    - `neout::Int`: Number of outboard engines.
    - `dyeout::Float64`: Distance between engines and the wingtip.
    - `neinn::Int`: Number of inboard engines.
    - `dyeinn::Float64`: Distance between engines and the wing root.
    - `Winn::Float64`: Weight of inboard engines.
    - `Wout::Float64`: Weight of outboard engines.
    - `dyWinn::Float64`: Weight distribution of inboard engines.
    - `dyWout::Float64`: Weight distribution of outboard engines.
    - `sweep::Float64`: Sweep angle in degrees.
    - `wbox::Float64`: Width of the wing box.
    - `tail.outboard.chord_thickness::Float64`: Height of the wing box at the root.
    - `tail.outboard.chord_thickness::Float64`: Height of the wing box at the strut attach point.
    - `rh::Float64`: Fractional height of the wing box.
    - `fLt::Float64`: Factor applied to the tip load.
    - `tauweb::Float64`: Web material shear strength.
    - `sigcap::Float64`: Cap material axial compressive strength.
    - `sigstrut::Float64`: Strut material axial compressive strength.
    - `Ecap::Float64`: Cap material Young's modulus.
    - `Eweb::Float64`: Web material Young's modulus.
    - `Gcap::Float64`: Cap material shear modulus.
    - `Gweb::Float64`: Web material shear modulus.
    - `rhoweb::Float64`: Density of the web material.
    - `rhocap::Float64`: Density of the cap material.
    - `rhostrut::Float64`: Density of the strut material.
    - `rhofuel::Float64`: Density of the fuel.
    
    **Outputs:**
    - `Ss::Float64`: Outboard section shear load.
    - `Ms::Float64`: Outboard section moment.
    - `tbwebs::Float64`: Web thickness at the strut attach point.
    - `tbcaps::Float64`: Cap thickness at the strut attach point.
    - `EIcs::Float64`: Combined cap and web bending stiffness at the strut attach point.
    - `EIns::Float64`: Combined cap and web normal stiffness at the strut attach point.
    - `GJs::Float64`: Combined cap and web shear stiffness at the strut attach point.
    - `So::Float64`: Inboard section shear load.
    - `Mo::Float64`: Inboard section moment.
    - `tbwebo::Float64`: Web thickness at the wing root.
    - `tbcapo::Float64`: Cap thickness at the wing root.
    - `EIco::Float64`: Combined cap and web bending stiffness at the wing root.
    - `EIno::Float64`: Combined cap and web normal stiffness at the wing root.
    - `GJo::Float64`: Combined cap and web shear stiffness at the wing root.
    - `Astrut::Float64`: Strut axial force.
    - `lsp::Float64`: Strut length.
    - `cosLs::Float64`: Cosine of the sweep angle at the strut attach point.
    - `Wscen::Float64`: Weight of center section (inboard of the strut).
    - `Wsinn::Float64`: Weight of the inner section.
    - `Wsout::Float64`: Weight of the outer section.
    - `dxWsinn::Float64`: Lateral distribution of inner section weight.
    - `dxWsout::Float64`: Lateral distribution of outer section weight.
    - `dyWsinn::Float64`: Vertical distribution of inner section weight.
    - `dyWsout::Float64`: Vertical distribution of outer section weight.
    - `Wfcen::Float64`: Weight of center section fuel.
    - `Wfinn::Float64`: Weight of the inner section fuel.
    - `Wfout::Float64`: Weight of the outer section fuel.
    - `dxWfinn::Float64`: Lateral distribution of inner section fuel weight.
    - `dxWfout::Float64`: Lateral distribution of outer section fuel weight.
    - `dyWfinn::Float64`: Vertical distribution of inner section fuel weight.
    - `dyWfout::Float64`: Vertical distribution of outer section fuel weight.
    - `Wweb::Float64`: Weight of the wing web.
    - `Wcap::Float64`: Weight of the wing cap.
    - `Wstrut::Float64`: Weight of the strut.
    - `dxWweb::Float64`: Lateral distribution of web weight.
    - `dxWcap::Float64`: Lateral distribution of cap weight.
    - `dxWstrut::Float64`: Lateral distribution of strut weight.

See [Geometry](@ref geometry),  [Wing/Tail Structures](@ref wingtail), and Section 2.7  of the [TASOPT Technical Description](@ref dreladocs). 
"""

function surft!(tail, po, lambdas,gammat,gammas,fLt,
            tauweb,sigcap,Ecap,rhoweb,rhocap, b = tail.outboard.b)

       Eweb = Ecap
       Gcap = Ecap * 0.5 / (1.0 + 0.3)
       Gweb = Ecap * 0.5 / (1.0 + 0.3)

       cosL = cos(deg2rad(tail.layout.sweep))
       sinL = sin(deg2rad(tail.layout.sweep))

      # Calculate non-dim span coordinate at span break and root (ηs and ηo resp.)
      etao = tail.layout.box_halfspan/b
      etas = tail.layout.box_halfspan/b

      cop = tail.layout.chord*cosL
      # Tip roll off Lift (modeled as a point load) and it's moment about ηs
      dLt = fLt*po*tail.layout.chord*gammat*tail.outboard.λ
      dMt = dLt*0.5*b*(1.0-etas)

      havgo = tail.outboard.chord_thickness * (1.0 - (1.0-tail.layout.hweb_to_hbox)/3.0)

      hrmso = tail.outboard.chord_thickness * sqrt(1.0 - (1.0-tail.layout.hweb_to_hbox)/1.5 + (1.0-tail.layout.hweb_to_hbox)^2 / 5.0)
      hrmss = tail.outboard.chord_thickness * sqrt(1.0 - (1.0-tail.layout.hweb_to_hbox)/1.5 + (1.0-tail.layout.hweb_to_hbox)^2 / 5.0)

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
      tbwebs = Ss*0.5/(cs^2*tauweb*tail.layout.hweb_to_hbox*tail.outboard.chord_thickness*cosL^2)
      con  = Ms*6.0*tail.outboard.chord_thickness/(cs^3*sigcap*tail.layout.box_width*cosL^4)
      tbcaps = 0.5*(hrmss - (hrmss^3-con)^(1.0/3.0))
      Abcaps = 2.0*tbcaps*tail.layout.box_width
      Abwebs = 2.0*tbwebs*tail.layout.hweb_to_hbox*tail.outboard.chord_thickness

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
    tbwebo = So*0.5/(tail.layout.chord^2*tauweb*tail.layout.hweb_to_hbox*tail.outboard.chord_thickness*cosL^2)
    con  = Mo*6.0*tail.outboard.chord_thickness/(tail.layout.chord^3*sigcap*tail.layout.box_width*cosL^4)
    tbcapo = 0.5*(hrmso - (hrmso^3-con)^(1.0/3.0))
    Abcapo = 2.0*tbcapo*tail.layout.box_width
    Abwebo = 2.0*tbwebo*tail.layout.hweb_to_hbox*tail.outboard.chord_thickness


      EIco = Ecap*cop^4 * (hrmso^3 - (hrmso-2.0*tbcapo)^3)*tail.layout.box_width/12.0 +
	 Eweb*cop^4 * tbwebo*(tail.layout.hweb_to_hbox*tail.outboard.chord_thickness)^3 / 6.0
      EIno = Ecap*cop^4 * tbcapo*    tail.layout.box_width  ^3 / 6.0 +
	 Eweb*cop^4 * tbwebo*tail.layout.hweb_to_hbox*tail.outboard.chord_thickness * 0.5*tail.layout.box_width^2
      GJo = cop^4 * 2.0*((tail.layout.box_width-tbwebo)*(havgo-tbcapo))^2 /
	 (  (tail.layout.hweb_to_hbox*tail.outboard.chord_thickness-tbcapo)/(Gweb*tbwebo) +
	 (   tail.layout.box_width -tbwebo)/(Gcap*tbcapo) )


      Vcen = tail.layout.chord^2*b *  etao / 2.0

      Vinn = tail.layout.chord^2*b * (etas-etao) *
	 (1.0 + lambdas + lambdas^2)/6.0 *
	 cosL
      Vout = tail.layout.chord^2*b * (1.0 -etas) *
	 (lambdas^2 + lambdas*tail.outboard.λ + tail.outboard.λ^2)/6.0 *
	 cosL

      dxVinn = tail.layout.chord^2*b^2 * (etas-etao)^2 *
	 (1.0 + 2.0*lambdas + 3.0*lambdas^2)/48.0 *
	 sinL
      dxVout = tail.layout.chord^2*b^2 * (1.0 -etas)^2 *
	 (lambdas^2 + 2.0*lambdas*tail.outboard.λ + 3.0*tail.outboard.λ^2)/48.0 *
	 sinL +
	  tail.layout.chord^2*b^2 * (etas-etao)*(1.0 -etas) *
	 (lambdas^2 + lambdas*tail.outboard.λ + tail.outboard.λ^2)/12.0 *
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
	 (lambdas^2 + lambdas*tail.outboard.λ + tail.outboard.λ^2)/6.0 *
	 cosL

        Whtail = tail.ntails * (Wscen + Wsinn + Wsout) * (1.0 + tail.weight_fraction_added)
        dxWhtail = tail.ntails * (dxWsinn + dxWsout) * (1.0 + tail.weight_fraction_added)
        tail.weight = Whtail
        tail.dxW = dxWhtail

     tail.thickness_web, tail.thickness_cap, tail.EI_bending, tail.EI_normal, tail.GJ =  tbwebo,tbcapo,EIco,EIno,GJo
end # surft!



