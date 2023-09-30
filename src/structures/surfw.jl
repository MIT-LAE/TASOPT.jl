"""
       surfw(gee, po, b, bs, bo, co, zs,
              lambdat, lambdas, gammat, gammas,
              Nload, iwplan, We, neout, dyeout, neinn, dyeinn,
              Winn, Wout, dyWinn, dyWout,
              sweep, wbox, hboxo, hboxs, rh, fLt,
              tauweb, sigcap, sigstrut, Ecap, Eweb, Gcap, Gweb,
              rhoweb, rhocap, rhostrut, rhofuel)

Calculates Wing or Tail loads, stresses, weights of individual wing sections.
Also returns the material gauges, torsional and bending stiffness.

Inputs :  

       gee ... rhostrut

Outputs:  

       Ss, Ms, tbwebs, tbcaps, EIcs, EIns, GJs,
              So, Mo, tbwebo, tbcapo, EIco, EIno, GJo,
              Astrut, lsp, cosLs,
              Wscen, Wsinn, Wsout, dxWsinn, dxWsout, dyWsinn, dyWsout,
              Wfcen, Wfinn, Wfout, dxWfinn, dxWfout, dyWfinn, dyWfout,
              Wweb, Wcap, Wstrut,
              dxWweb, dxWcap, dxWstrut

"""
function surfw(gee,po,b,bs,bo,co,zs,
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

      # Inboard engines


#---- strut-attach shear,moment from outer-wing loading. Note added term to account for outboard engines
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

#      if(iwplan==0) 
#c----- no strut, no engine
#       ls = 0.
#       Tstrut = 0.
#       Rstrut = 0.
#       Pstrut = 0.
#       So = Ss
#       Mo = Ms
#       tbwebo = tbwebs
#       tbcapo = tbcaps
#       Abwebo = Abwebs
#       Abcapo = Abcaps
#       lsp = 0.
#       Tstrutp = 0.
#       cosLs = 1.0

      if(iwplan==0 || iwplan==1) 
#----- no strut, with or without engine at etas
       ls = 0.
       Tstrut = 0.
       Rstrut = 0.
       Pstrut = 0.
# Original from TASOPT
       # So = Ss - Nload*We +
	#  0.25*po*b*(1.0+gammas)*(etas-etao) -
	#  Nload*Winn
       # Mo = Ms + (Ss-Nload*We)*0.5*b*(etas-etao) +
	#  (1.0/24.0)*po*b^2*(1.0+2.0*gammas)*(etas-etao)^2 -
	#  Nload*dyWinn
# Modifed to account for bending relief from multiple engines:
       So = Ss - Nload*neinn*We +
	 0.25*po*b*(1.0+gammas)*(etas-etao) -
	 Nload*Winn
       Mo = Ms + Ss*0.5*b*(etas-etao) +
	 (1.0/24.0)*po*b^2*(1.0+2.0*gammas)*(etas-etao)^2 -
	 Nload*dyWinn - Nload*neinn*We*dyeinn

#----- limit So,Mo to Ss,Ms, which might be needed with heavy outboard engine
#-      (rules out negatively-tapered structure, deemed not feasible for downloads)
       So = max( So , Ss )
       Mo = max( Mo , Ms )

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


