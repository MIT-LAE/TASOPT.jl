"""
    surfw(po, b, bs, bo, co, zs,
        lambdat, lambdas, gammat, gammas,
        Nload, iwplan, We, neout, dyeout, neinn, dyeinn,
        Winn, Wout, dyWinn, dyWout,
        sweep, wbox, hboxo, hboxs, rh, fLt,
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
    - `lambdat::Float64`: Tip chord ratio (tip chord / root chord).
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
    - `hboxo::Float64`: Height of the wing box at the root.
    - `hboxs::Float64`: Height of the wing box at the strut attach point.
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

See [here](@ref wingtail) or Section 2.7  of TASOPT docs. 
"""
function surfw(po,b,bs,bo,co,zs,
	lambdat,lambdas,gammat,gammas,
	Nload,iwplan,We,neout, dyeout, neinn, dyeinn,
	Winn,Wout,dyWinn,dyWout,
	sweep,wbox,hboxo,hboxs,rh, fLt,
	tauweb,sigcap,sigstrut,Ecap,Eweb,Gcap,Gweb,
	rhoweb,rhocap,rhostrut,rhofuel)

      cosL = cos(sweep*pi/180)
      sinL = sin(sweep*pi/180)

      # Calculate non-dim span coordinate at span break and root (ηs and ηo resp.)
      etao = bo/b
      etas = bs/b

      cop = co*cosL
      csp = co*cosL*lambdas
      # Tip roll off Lift (modeled as a point load) and it's moment about ηs
      dLt = fLt*po*co*gammat*lambdat
      dMt = dLt*0.5*b*(1.0-etas)

      havgo = hboxo * (1.0 - (1.0-rh)/3.0)
      havgs = hboxs * (1.0 - (1.0-rh)/3.0)

      hrmso = hboxo * sqrt(1.0 - (1.0-rh)/1.5 + (1.0-rh)^2 / 5.0)
      hrmss = hboxs * sqrt(1.0 - (1.0-rh)/1.5 + (1.0-rh)^2 / 5.0)

# Outboard section:
#---- strut-attach shear,moment from outer-wing loading. 
#     Note added term to account for any outboard engines.
#     If neout = 0 this simplifies to Drela's version which assumes engine
#     fixed at ηs locations only.
      Ss = (po*b   / 4.0)*(gammas+    gammat)*(1.0-etas) +
	 dLt - Nload*Wout - Nload*neout*We
      Ms = (po*b^2/24.0)*(gammas+2.0*gammat)*(1.0-etas)^2 +
	 dMt - Nload*dyWout - Nload*neout*We*dyeout

#---- size strut-attach station at etas
      cs = co*lambdas
      tbwebs = Ss*0.5/(cs^2*tauweb*rh*hboxs*cosL^2)
      con  = Ms*6.0*hboxs/(cs^3*sigcap*wbox*cosL^4)
      tbcaps = 0.5*(hrmss - (hrmss^3-con)^(1.0/3.0))
      Abcaps = 2.0*tbcaps*wbox
      Abwebs = 2.0*tbwebs*rh*hboxs

      EIcs = Ecap*csp^4 * (hrmss^3 - (hrmss-2.0*tbcaps)^3)*wbox/12.0 +
	 Eweb*csp^4 * tbwebs*(rh*hboxs)^3 / 6.0
      EIns = Ecap*csp^4 * tbcaps*    wbox  ^3 / 6.0 +
	 Eweb*csp^4 * tbwebs*rh*hboxs * 0.5*wbox^2
      GJs = csp^4 * 2.0*((wbox-tbwebs)*(havgs-tbcaps))^2 /
	 (  (rh*hboxs-tbcaps)/(Gweb*tbwebs) +
	 (   wbox -tbwebs)/(Gcap*tbcaps) )

# Inboard Section:
      if(iwplan==0 || iwplan==1) 
#----- no strut, with or without engine at etas
       ls = 0.
       Tstrut = 0.
       Rstrut = 0.
       Pstrut = 0.

# Modifed to account for bending relief from multiple engines.
# dyeinn allows engine to be in locations other than ηs
       So = Ss - Nload*neinn*We +
	 0.25*po*b*(1.0+gammas)*(etas-etao) -
	 Nload*Winn
       Mo = Ms + Ss*0.5*b*(etas-etao) +
	 (1.0/24.0)*po*b^2*(1.0+2.0*gammas)*(etas-etao)^2 -
	 Nload*dyWinn - Nload*neinn*We*dyeinn

#----- limit So,Mo to Ss,Ms, which might be needed with heavy outboard engine
#-      (rules out negatively-tapered structure, deemed not feasible for downloads)
       So = max(So ,Ss)
       Mo = max(Mo ,Ms)

#----- size root station at etao
       tbwebo = So*0.5/(co^2*tauweb*rh*hboxo*cosL^2)
       con  = Mo*6.0*hboxo/(co^3*sigcap*wbox*cosL^4)
       tbcapo = 0.5*(hrmso - (hrmso^3-con)^(1.0/3.0))
       Abcapo = 2.0*tbcapo*wbox
       Abwebo = 2.0*tbwebo*rh*hboxo

       lsp = 0.
       Tstrutp = 0.
       cosLs = 1.0

      else
#----- strut present
       ls = sqrt(zs^2 + (0.5*b*(etas-etao))^2)
       Rstrut = (po*b/12.0)*(etas-etao)*(1.0+2.0*gammas) + Ss
       Tstrut = Rstrut*ls/zs
#c     Pstrut = Rstrut*0.5*b*(etas-etao)/zs

#----- inboard shear,moment used for sparbox sizing
       So = Ss
       Mo = Ms
#
#----- size inboard station at etao
       tbwebo = So*0.5/(co^2*tauweb*rh*hboxo*cosL^2)
       con  = Mo*6.0*hboxo/(co^3*sigcap*wbox*cosL^4)
       tbcapo = 0.5*(hrmso - (hrmso^3-con)^(1.0/3.0))
       Abcapo = 2.0*tbcapo*wbox
       Abwebo = 2.0*tbwebo*rh*hboxo

#----- total strut length, tension
       lsp = sqrt(zs^2 + (0.5*b*(etas-etao)/cosL)^2)
       Tstrutp = Tstrut*lsp/ls
       cosLs = ls/lsp
      
      end

      EIco = Ecap*cop^4 * (hrmso^3 - (hrmso-2.0*tbcapo)^3)*wbox/12.0 +
	 Eweb*cop^4 * tbwebo*(rh*hboxo)^3 / 6.0
      EIno = Ecap*cop^4 * tbcapo*    wbox  ^3 / 6.0 +
	 Eweb*cop^4 * tbwebo*rh*hboxo * 0.5*wbox^2
      GJo = cop^4 * 2.0*((wbox-tbwebo)*(havgo-tbcapo))^2 /
	 (  (rh*hboxo-tbcapo)/(Gweb*tbwebo) +
	 (   wbox -tbwebo)/(Gcap*tbcapo) )


      Abfuels = (wbox-2.0*tbwebs)*(havgs-2.0*tbcaps)
      Abfuelo = (wbox-2.0*tbwebo)*(havgo-2.0*tbcapo)

      Astrut = Tstrutp/sigstrut

      Vcen = co^2*b *  etao / 2.0

      Vinn = co^2*b * (etas-etao) *
	 (1.0 + lambdas + lambdas^2)/6.0 *
	 cosL
      Vout = co^2*b * (1.0 -etas) *
	 (lambdas^2 + lambdas*lambdat + lambdat^2)/6.0 *
	 cosL

      dxVinn = co^2*b^2 * (etas-etao)^2 *
	 (1.0 + 2.0*lambdas + 3.0*lambdas^2)/48.0 *
	 sinL
      dxVout = co^2*b^2 * (1.0 -etas)^2 *
	 (lambdas^2 + 2.0*lambdas*lambdat + 3.0*lambdat^2)/48.0 *
	 sinL +
	  co^2*b^2 * (etas-etao)*(1.0 -etas) *
	 (lambdas^2 + lambdas*lambdat + lambdat^2)/12.0 *
	 sinL

      dyVinn = co^2*b^2 * (etas-etao)^2 *
	 (1.0 + 2.0*lambdas + 3.0*lambdas^2)/48.0 *
	 cosL
      dyVout = co^2*b^2 * (1.0 -etas)^2 *
	 (lambdas^2 + 2.0*lambdas*lambdat + 3.0*lambdat^2)/48.0 *
	 cosL

#---- set chord^2 weighted average areas for inner panel
      Abcapi = (Abcapo + Abcaps*lambdas^2)/(1.0+lambdas^2)
      Abwebi = (Abwebo + Abwebs*lambdas^2)/(1.0+lambdas^2)

      Wscen   = (rhocap*Abcapo + rhoweb*Abwebo)*gee*Vcen
      Wsinn   = (rhocap*Abcapi + rhoweb*Abwebi)*gee*Vinn
      Wsout   = (rhocap*Abcaps + rhoweb*Abwebs)*gee*Vout

      dxWsinn = (rhocap*Abcapi + rhoweb*Abwebi)*gee*dxVinn
      dxWsout = (rhocap*Abcaps + rhoweb*Abwebs)*gee*dxVout

      dyWsinn = (rhocap*Abcapi + rhoweb*Abwebi)*gee*dyVinn
      dyWsout = (rhocap*Abcaps + rhoweb*Abwebs)*gee*dyVout


      Abfueli= (Abfuelo + Abfuels*lambdas^2)/(1.0+lambdas^2)

      Wfcen   = rhofuel*Abfuelo *gee*Vcen
      Wfinn   = rhofuel*Abfueli *gee*Vinn
      Wfout   = rhofuel*Abfuels *gee*Vout

      dxWfinn = rhofuel*Abfueli *gee*dxVinn
      dxWfout = rhofuel*Abfuels *gee*dxVout

      dyWfinn = rhofuel*Abfueli *gee*dyVinn
      dyWfout = rhofuel*Abfuels *gee*dyVout

      Wcap   = 2.0*rhocap*gee*(Abcapo*Vcen + Abcapi*Vinn + Abcaps*Vout)
      Wweb   = 2.0*rhoweb*gee*(Abwebo*Vcen + Abwebi*Vinn + Abwebs*Vout)

      dxWcap = 2.0*rhocap*gee*( Abcapi*dxVinn + Abcaps*dxVout )
      dxWweb = 2.0*rhoweb*gee*( Abwebi*dxVinn + Abwebs*dxVout )

      Wstrut   = 2.0*rhostrut*gee*Astrut*lsp
      dxWstrut = Wstrut * 0.25*b*(etas-etao) * sinL/cosL

      Vout = co^2*b * (1.0 -etas) *
	 (lambdas^2 + lambdas*lambdat + lambdat^2)/6.0 *
	 cosL


      return	Ss,Ms,tbwebs,tbcaps,EIcs,EIns,GJs,
	So,Mo,tbwebo,tbcapo,EIco,EIno,GJo,
	Astrut,lsp,cosLs,
	Wscen,Wsinn,Wsout,dxWsinn,dxWsout,dyWsinn,dyWsout,
	Wfcen,Wfinn,Wfout,dxWfinn,dxWfout,dyWfinn,dyWfout,
	Wweb,  Wcap,  Wstrut,
	dxWweb,dxWcap,dxWstrut

end # surfw


