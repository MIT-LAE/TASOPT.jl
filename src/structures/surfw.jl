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

tauweb,sigcap,sigstrut = wing.inboard.webs.material.τmax * sigfac, wing.inboard.caps.σ * sigfac, wing.strut.material.σmax * sigfac

Eweb = wing.inboard.caps.material.E
Gcap = wing.inboard.caps.material.E * 0.5 / (1.0 + 0.3)
Gweb = wing.inboard.caps.material.E * 0.5 / (1.0 + 0.3)

cosL = cos(wing.layout.sweep*pi/180)
sinL = sin(wing.layout.sweep*pi/180)

# Calculate non-dim span coordinate at span break and root (ηs and ηo resp.)
etao = wing.layout.box_halfspan/wing.outboard.layout.b
etas = wing.inboard.layout.b/wing.outboard.layout.b

cop = wing.layout.chord*cosL
csp = wing.layout.chord*cosL*wing.inboard.layout.λ
# Tip roll off Lift (modeled as a point load) and it's moment about ηs
dLt = fLt*po*wing.layout.chord*gammat*wing.outboard.layout.λ
dMt = dLt*0.5*wing.outboard.layout.b*(1.0-etas)

havgo = wing.inboard.layout.chord_thickness * (1.0 - (1.0-wing.layout.hweb_to_hbox)/3.0)
havgs = wing.outboard.layout.chord_thickness * (1.0 - (1.0-wing.layout.hweb_to_hbox)/3.0)

hrmso = wing.inboard.layout.chord_thickness * sqrt(1.0 - (1.0-wing.layout.hweb_to_hbox)/1.5 + (1.0-wing.layout.hweb_to_hbox)^2 / 5.0)
hrmss = wing.outboard.layout.chord_thickness * sqrt(1.0 - (1.0-wing.layout.hweb_to_hbox)/1.5 + (1.0-wing.layout.hweb_to_hbox)^2 / 5.0)

# Outboard section:
#---- strut-attach shear,moment from outer-wing loading. 
#     Note added term to account for any outboard engines.
#     If neout = 0 this simplifies to Drela's version which assumes engine
#     fixed at ηs locations only.
wing.outboard.max_shear_load = (po*wing.outboard.layout.b   / 4.0)*(gammas+    gammat)*(1.0-etas) +
 dLt - Nload*wing.outboard.weight - Nload*neout*We
 wing.outboard.moment = (po*wing.outboard.layout.b^2/24.0)*(gammas+2.0*gammat)*(1.0-etas)^2 +
 dMt - Nload*wing.outboard.dyW - Nload*neout*We*dyeout

#---- size strut-attach station at etas
cs = wing.layout.chord*wing.inboard.layout.λ
tbwebs = wing.outboard.max_shear_load*0.5/(cs^2*tauweb*wing.layout.hweb_to_hbox*wing.outboard.layout.chord_thickness*cosL^2)
con  = wing.outboard.moment*6.0*wing.outboard.layout.chord_thickness/(cs^3*sigcap*wing.layout.box_width*cosL^4)
tbcaps = 0.5*(hrmss - (hrmss^3-con)^(1.0/3.0))
Abcaps = 2.0*tbcaps*wing.layout.box_width
Abwebs = 2.0*tbwebs*wing.layout.hweb_to_hbox*wing.outboard.layout.chord_thickness

wing.outboard.web_cap.EI_bending = wing.inboard.caps.material.E*csp^4 * (hrmss^3 - (hrmss-2.0*tbcaps)^3)*wing.layout.box_width/12.0 +
 Eweb*csp^4 * tbwebs*(wing.layout.hweb_to_hbox*wing.outboard.layout.chord_thickness)^3 / 6.0
 wing.outboard.web_cap.EI_normal = wing.inboard.caps.material.E*csp^4 * tbcaps*    wing.layout.box_width  ^3 / 6.0 +
 Eweb*csp^4 * tbwebs*wing.layout.hweb_to_hbox*wing.outboard.layout.chord_thickness * 0.5*wing.layout.box_width^2
wing.outboard.web_cap.GJ = csp^4 * 2.0*((wing.layout.box_width-tbwebs)*(havgs-tbcaps))^2 /
 (  (wing.layout.hweb_to_hbox*wing.outboard.layout.chord_thickness-tbcaps)/(Gweb*tbwebs) +
 (   wing.layout.box_width -tbwebs)/(Gcap*tbcaps) )

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
 0.25*po*wing.outboard.layout.b*(1.0+gammas)*(etas-etao) -
 Nload*wing.inboard.weight
Mo = wing.outboard.moment + wing.outboard.max_shear_load*0.5*wing.outboard.layout.b*(etas-etao) +
 (1.0/24.0)*po*wing.outboard.layout.b^2*(1.0+2.0*gammas)*(etas-etao)^2 -
 Nload*wing.inboard.dyW - Nload*neinn*We*dyeinn

#----- limit So,Mo to Ss,Ms, which might be needed with heavy outboard engine
#-      (rules out negatively-tapered structure, deemed not feasible for downloads)
wing.inboard.max_shear_load = max(So ,wing.outboard.max_shear_load)
wing.inboard.moment = max(Mo ,wing.outboard.moment)

#----- size root station at etao
tbwebo = So*0.5/(wing.layout.chord^2*tauweb*wing.layout.hweb_to_hbox*wing.inboard.layout.chord_thickness*cosL^2)
con  = Mo*6.0*wing.inboard.layout.chord_thickness/(wing.layout.chord^3*sigcap*wing.layout.box_width*cosL^4)
tbcapo = 0.5*(hrmso - (hrmso^3-con)^(1.0/3.0))
Abcapo = 2.0*tbcapo*wing.layout.box_width
Abwebo = 2.0*tbwebo*wing.layout.hweb_to_hbox*wing.inboard.layout.chord_thickness

lsp = 0.
Tstrutp = 0.
wing.strut.cos_lambda = 1.0

else
#----- strut present
ls = sqrt(wing.strut.z^2 + (0.5*wing.outboard.layout.b*(etas-etao))^2)
Rstrut = (po*wing.outboard.layout.b/12.0)*(etas-etao)*(1.0+2.0*gammas) + wing.outboard.max_shear_load
Tstrut = Rstrut*ls/wing.strut.z
#c     Pstrut = Rstrut*0.5*wing.outboard.layout.b*(etas-etao)/zs

#----- inboard shear,moment used for sparbox sizing
wing.inboard.max_shear_load = wing.outboard.max_shear_load
wing.inboard.moment = wing.outboard.moment
#
#----- size inboard station at etao
tbwebo = wing.inboard.max_shear_load*0.5/(wing.layout.chord^2*tauweb*wing.layout.hweb_to_hbox*wing.inboard.layout.chord_thickness*cosL^2)
con  = wing.inboard.moment*6.0*wing.inboard.layout.chord_thickness/(wing.layout.chord^3*sigcap*wing.layout.box_width*cosL^4)
tbcapo = 0.5*(hrmso - (hrmso^3-con)^(1.0/3.0))
Abcapo = 2.0*tbcapo*wing.layout.box_width
Abwebo = 2.0*tbwebo*wing.layout.hweb_to_hbox*wing.inboard.layout.chord_thickness

#----- total strut length, tension
lsp = sqrt(wing.strut.z^2 + (0.5*wing.outboard.layout.b*(etas-etao)/cosL)^2)
Tstrutp = Tstrut*lsp/ls
wing.strut.cos_lambda = ls/lsp

end

wing.inboard.web_cap.EI_bending = wing.inboard.caps.material.E*cop^4 * (hrmso^3 - (hrmso-2.0*tbcapo)^3)*wing.layout.box_width/12.0 +
 Eweb*cop^4 * tbwebo*(wing.layout.hweb_to_hbox*wing.inboard.layout.chord_thickness)^3 / 6.0
 wing.inboard.web_cap.EI_normal = wing.inboard.caps.material.E*cop^4 * tbcapo*    wing.layout.box_width  ^3 / 6.0 +
 Eweb*cop^4 * tbwebo*wing.layout.hweb_to_hbox*wing.inboard.layout.chord_thickness * 0.5*wing.layout.box_width^2
wing.inboard.web_cap.GJ = cop^4 * 2.0*((wing.layout.box_width-tbwebo)*(havgo-tbcapo))^2 /
 (  (wing.layout.hweb_to_hbox*wing.inboard.layout.chord_thickness-tbcapo)/(Gweb*tbwebo) +
 (   wing.layout.box_width -tbwebo)/(Gcap*tbcapo) )


Abfuels = (wing.layout.box_width-2.0*tbwebs)*(havgs-2.0*tbcaps)
Abfuelo = (wing.layout.box_width-2.0*tbwebo)*(havgo-2.0*tbcapo)

wing.strut.axial_force = Tstrutp/sigstrut

Vcen = wing.layout.chord^2*wing.outboard.layout.b*  etao / 2.0

Vinn = wing.layout.chord^2*wing.outboard.layout.b* (etas-etao) *
 (1.0 + wing.inboard.layout.λ + wing.inboard.layout.λ^2)/6.0 *
 cosL
Vout = wing.layout.chord^2*wing.outboard.layout.b* (1.0 -etas) *
 (wing.inboard.layout.λ^2 + wing.inboard.layout.λ*wing.outboard.layout.λ + wing.outboard.layout.λ^2)/6.0 *
 cosL

dxVinn = wing.layout.chord^2*wing.outboard.layout.b^2 * (etas-etao)^2 *
 (1.0 + 2.0*wing.inboard.layout.λ + 3.0*wing.inboard.layout.λ^2)/48.0 *
 sinL
dxVout = wing.layout.chord^2*wing.outboard.layout.b^2 * (1.0 -etas)^2 *
 (wing.inboard.layout.λ^2 + 2.0*wing.inboard.layout.λ*wing.outboard.layout.λ + 3.0*wing.outboard.layout.λ^2)/48.0 *
 sinL +
  wing.layout.chord^2*wing.outboard.layout.b^2 * (etas-etao)*(1.0 -etas) *
 (wing.inboard.layout.λ^2 + wing.inboard.layout.λ*wing.outboard.layout.λ + wing.outboard.layout.λ^2)/12.0 *
 sinL

dyVinn = wing.layout.chord^2*wing.outboard.layout.b^2 * (etas-etao)^2 *
 (1.0 + 2.0*wing.inboard.layout.λ + 3.0*wing.inboard.layout.λ^2)/48.0 *
 cosL
dyVout = wing.layout.chord^2*wing.outboard.layout.b^2 * (1.0 -etas)^2 *
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
wing.strut.dxW = wing.strut.weight * 0.25*wing.outboard.layout.b*(etas-etao) * sinL/cosL

Vout = wing.layout.chord^2*wing.outboard.layout.b* (1.0 -etas) *
 (wing.inboard.layout.λ^2 + wing.inboard.layout.λ*wing.outboard.layout.λ + wing.outboard.layout.λ^2)/6.0 *
 cosL

 wing.inboard.caps.thickness = tbcapo
 wing.inboard.webs.thickness = tbwebo
 wing.outboard.caps.thickness = tbcaps
 wing.outboard.webs.thickness = tbwebs

 fwadd = wing_additional_weight(wing)
 Wwing = 2.0 * (Wscen + Wsinn + Wsout) * (1.0 + fwadd)
 wing.dxW = 2.0 * (dxWsinn + dxWsout) * (1.0 + fwadd)


#  wing.outboard.max_shear_load, wing.outboard.moment, wing.outboard.webs.thickness, wing.outboard.caps.thickness, wing.outboard.web_cap.EI_bending, wing.outboard.web_cap.EI_normal, wing.outboard.web_cap.GJ,
#  wing.inboard.max_shear_load, wing.inboard.moment, wing.inboard.webs.thickness, wing.inboard.caps.thickness, wing.inboard.web_cap.EI_bending, wing.inboard.web_cap.EI_normal, wing.inboard.web_cap.GJ,
#  wing.strut.axial_force, lstrutp, wing.strut.cos_lambda,
#  Wscen, Wsinn, Wsout, dxWsinn, dxWsout, dyWsinn, dyWsout,
#  Wfcen, Wfinn, Wfout, dxWfinn, dxWfout, dyWfinn, dyWfout,
#  wing.inboard.webs.weight.W, wing.inboard.caps.weight.W, wing.strut.weight,
#  dxWweb, dxWcap, wing.strut.dxW

#       return	Ss,Ms,tbwebs,tbcaps,EIcs,EIns,GJs,
# 	So,Mo,tbwebo,tbcapo,EIco,EIno,GJo,




return Wwing,Wsinn,Wsout,dyWsinn,dyWsout,Wfcen,Wfinn,Wfout,dxWfinn,dxWfout,dyWfinn,dyWfout,lsp



end # surfw


