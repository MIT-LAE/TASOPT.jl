
"""
      balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; Ldebug)

Computes the aircraft's center of gravity (`xCG`), center of pressure (`xCP`), and neutral point (`xNP`) based on payload, fuel distribution, and trim adjustments.
Makes one of three (or none) changes to achieve pitch trim. Formerly, `balance()`.

**Description**
This routine performs a CG and stability analysis for a given aircraft configuration. It calculates the **total weight and moment** by accounting for:
- Payload distribution (`rpay`, `ξpay`).
- Fuel distribution (`rfuel`).
- Structural components (fuselage, wing, tail, landing gear).
- Trim adjustments (`opt_trim_var`), which modify horizontal tail lift, area, or wing box location.

The routine computes the **neutral point (`xNP`), indicating the aircraft's longitudinal static stability, and may achieve pitch trim by adjusting for one of the following:

| `opt_trim_var` | Adjustment Method |
|----------------|-------------------|
| `"none"`       | No adjustments. Only calculates and returns the neutral point (`xNP`) | 
| `"CL_htail"`   | Adjusts horizontal tail lift coefficient (`CLh`) |
| `"S_htail"`    | Adjusts horizontal tail area (`Sh`) |
| `"x_wingbox"`  | Adjusts wing box location (`xwbox`) |

!!! details "🔃 Inputs and Outputs"
      **Inputs**
      - `ac` : Aircraft object
      - `imission` : Mission index (used for unpacking mission-specific parameters; 1 is design mission).
      - `ip` : flight point index (used for aerodynamic/weight calculations).
      - `rfuel` : Fuel fraction.
      - `rpay` : Payload fraction.
      - `ξpay` : Payload distribution factor (0.0 = front-loaded, 1.0 = rear-loaded).
      - `opt_trim_var` : Variable to adjust to achieve pitch trim (`"none"` for no adjustments, `"CL_htail"` for htail lift coefficient, `"S_htail"` for htail area, `"x_wingbox"` for wing box location).
      - `Ldebug` : Optional debug flag (default: `false`). If `true`, prints debug information.

      **Outputs** 
      No explicit return values, but updates fields inside `para`. Namely:
      - `para[iaxCG]` : Computed center of gravity (`xCG`).
      - `para[iaxCP]` : Computed center of pressure (`xCP`).
      - `para[iaxNP]` : Computed neutral point (`xNP`).

**Notes**
- Uses [`CG_limits()`](@ref TASOPT.CG_limits) to compute CG limits (`xcgF`, `xcgB`).
- Uses [`cabin_centroid()`](@ref TASOPT.cabin_centroid) to determine cabin location.
- If there is fuel in the wings (`ac.options.has_wing_fuel`), it does not shift between CG cases.
- `xNP` is affected by engine placement (`xengcp`), aerodynamics (`CMw1`, `CMh1`), and fuel distribution.

"""
function balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; Ldebug::Bool = false)
      #Unpack aircraft
      parg, _, para, _, options, fuse, fuse_tank, wing, htail, vtail, _, landing_gear = unpack_ac(ac, imission, ip = ip)

      # Unpack weights
      Wpay = parg[igWpay]
      Wfuel = parg[igWfuel]
      Wfuse = fuse.weight
      Wwing = wing.weight
      Wstrut = wing.strut.weight
      Whtail = htail.weight
      Wvtail = vtail.weight
      Weng = parg[igWeng]

      xWfuel = parg[igxWfuel]

      Wtesys = parg[igWtesys]
      xWtesys = parg[igxWtesys]

      Wftank = parg[igWftank]
      xWftank = parg[igxWftank]

      nftanks = fuse_tank.tank_count #number of fuel tanks in fuselage
      lftank = parg[iglftank]

      # Use weight fractions to calculate weights of subsystems
      Whpesys = parg[igWMTO] * fuse.HPE_sys.W
      Wlgnose = landing_gear.nose_gear.weight.W
      Wlgmain = landing_gear.main_gear.weight.W

      xcabin,lcabin = cabin_centroid(nftanks,fuse,parg[igxftankaft],lftank)
      
      xpay = xcabin + (ξpay - 0.5) * lcabin * (1.0 - rpay)

      xwbox = wing.layout.box_x

      rfuelF, rfuelB, rpayF, rpayB, xcgF, xcgB = CG_limits(ac; Ldebug = Ldebug)

      #---- wing centroid offset from wingbox, assumed fixed in CG calculations
      dxwing = wing.layout.x - wing.layout.box_x

      #---- main LG offset from wingbox, assumed fixed in CG calculations
      dxlg = xcgB + landing_gear.main_gear.distance_CG_to_landing_gear - wing.layout.box_x

      S = wing.layout.S
      Sh = htail.layout.S
      co = wing.layout.root_chord
      coh = htail.layout.root_chord
      xhbox = htail.layout.box_x

      cma = wing.mean_aero_chord

      #---- for better convergence, 
      #      will assume that HT weight will scale with its area
      Sh1 = Sh

      #---- total weight, weight moment, and derivatives
      W = rpay * Wpay +
          rfuel * Wfuel +
          Wfuse +
          Wtesys +
          Wftank +
          Wwing +
          Wstrut +
          Whtail * Sh / Sh1 +
          Wvtail +
          Weng +
          Whpesys +
          Wlgnose +
          Wlgmain

      #  println("rpay = $rpay, rufel = $rfuel;
      #   rpay *Wpay  = $(rpay *Wpay ) 
      #  rfuel*Wfuel = $(rfuel*Wfuel ) 
      #  Wfuse  = $(Wfuse  ) 
      #  Wtesys = $(Wtesys ) 
      #  Wftank = $(Wftank ) 
      #  Wwing  = $(Wwing  ) 
      #  Wstrut = $(Wstrut ) 
      #  Whtail*Sh/Sh1 = $(Whtail*Sh/Sh1) 
      #  Wvtail  = $(Wvtail  ) 
      #  Weng    = $(Weng    ) 
      #  Whpesys = $(Whpesys ) 
      #  Wlgnose = $(Wlgnose ) 
      #  Wlgmain = $(Wlgmain )")

      W_Sh = Whtail / Sh1

      #Calculate fuel moment and derivatives wrt to wing box location 
      if !(options.has_wing_fuel) #If no wing fuel, therefore fuel is stored in the fuselage
            xWfuel = rfuel * xWfuel
            xWfuel_xwbox = 0.0
      else #if wing fuel
            xWfuel = rfuel * xWfuel
            xWfuel_xwbox = rfuel * Wfuel
      end

      xW = rpay * Wpay * xpay +
           xWfuel +
           fuse.moment + xWtesys + xWftank +
           Wwing * wing.layout.box_x + wing.dxW +
           Wstrut * wing.layout.box_x + wing.strut.dxW +
           (Whtail * htail.layout.box_x + htail.dxW) * Sh / Sh1 +
           Wvtail * vtail.layout.box_x + vtail.dxW +
           Weng * parg[igxeng] +
           Whpesys * fuse.HPE_sys.r.x +
           Wlgmain * (wing.layout.box_x + dxlg) + landing_gear.nose_gear.moment

      xW_xwbox = xWfuel_xwbox + Wwing + Wstrut + Wlgmain

      xW_Sh = (Whtail * htail.layout.box_x + htail.dxW) / Sh1

      #---- total aero moment and derivatives
      CMw0 = para[iaCMw0]
      CMw1 = para[iaCMw1]
      CMh0 = para[iaCMh0]
      CMh1 = para[iaCMh1]
      CL = para[iaCL]
      CLh = para[iaCLh]
      # println("before opt_trim_var $(para[iaCLh]), CL = $CL")
      CMVf1 = parg[igCMVf1]
      CLMf0 = parg[igCLMf0]

      cCM = co * CMw0 + (co * CMw1 - xwbox) * (CL - CLh * Sh / S) +
            coh * CMh0 * Sh / S + (coh * CMh1 - xhbox) * CLh * Sh / S +
            CMVf1 * (CL - CLMf0) / S

      cCM_xwbox = -(CL - CLh * Sh / S)
      cCM_Sh = (co * CMw1 - xwbox) * (-CLh / S) +
               coh * CMh0 / S + (coh * CMh1 - xhbox) * CLh / S
      cCM_CLh = (co * CMw1 - xwbox) * (-Sh / S) +
                (coh * CMh1 - xhbox) * Sh / S

      #---- current CP and CG location
      #      xCP = -cCM/CL
      #      xCG = xW/W

      #---- total-moment residual
      Res = cCM / CL + xW / W                 # Eq. 300
      Res_CLh = cCM_CLh / CL
      Res_Sh = cCM_Sh / CL + xW_Sh / W - (xW / W^2) * W_Sh
      Res_xwbox = cCM_xwbox / CL + xW_xwbox / W

      #      r[1] = Res
      #      a[1,1] = Res_CLh
      #      a[1,2] = Res_Sh
      #      a[1,3] = Res_xwbox

      #---- drive to Res=0 (pitch trim) by one of three changes...

      if compare_strings(opt_trim_var,"CL_htail")
            #----- ... adjust horizontal tail CLh
            delCLh = -Res / Res_CLh
            CLh = CLh + delCLh
            #  println("inside opt_trim_var = "CL_htail": $delCLh, $CLh")
            cCM = cCM + cCM_CLh * delCLh
            para[iaCLh] = CLh

      elseif compare_strings(opt_trim_var, "S_htail")
            #----- ... adjust horizontal tail area
            delSh = -Res / Res_Sh
            Sh = Sh + delSh
            cCM = cCM + cCM_Sh * delSh
            xW = xW + xW_Sh * delSh

            htail.layout.S = Sh

      elseif compare_strings(opt_trim_var, "x_wingbox")
            #----- ... adjust wing box location
            delxwbox = -Res / Res_xwbox
            xwbox = xwbox + delxwbox
            cCM = cCM + cCM_xwbox * delxwbox
            xW = xW + xW_xwbox * delxwbox

            wing.layout.box_x = xwbox
            wing.layout.x = xwbox + dxwing

      end

      #---- calculate neutral point
      dCLhdCL = parg[igdCLhdCL]

      dCLndCL = parg[igdCLndCL]
      neng = parg[igneng]
      xeng = parg[igxeng]
      dfan = parg[igdfan]
      Afan = 0.25 * pi * dfan^2
      xengcp = xeng - 0.25 * dfan

      xNP = (co * CMw1 - xwbox) * (dCLhdCL * Sh / S - 1.0) -
            (coh * CMh1 - xhbox) * dCLhdCL * Sh / S -
            CMVf1 / S -
            neng * xengcp * dCLndCL * Afan / S
      xNP_xwbox = -(dCLhdCL * Sh / S - 1.0)
      xNP_Sh = (co * CMw1 - xwbox) * dCLhdCL / S -
               (coh * CMh1 - xhbox) * dCLhdCL / S
      xNP_CLh = 0.0

      #---- x locations after one of the trim changes
      #-     (xCP,xCG will be equal if moments are linear)
      xCP = -cCM / CL
      xCG = xW / W

      #---- component pitching moments
      #      cCM  = co *CMw0      + (co *CMw1 - xwbox)*(CL - CLh*Sh/S) - xref*CL
      #     &     + coh*CMh0*Sh/S + (coh*CMh1 - xhbox)*      CLh*Sh/S
      #     &     + CMVf1*(CL-CLMf0)/S
      #
      #      cCM_CL  = (co*CMw1  - xwbox - xref)*(1.0 - CLh_CL*Sh/S)
      #     &        + (coh*CMh1 - xhbox - xref)*CLh_CL*Sh/S

      #     CMwing = (co *CMw0      + (co *CMw1-xwbox-xCG)*(CL-CLh*Sh/S))/cma
      #     CMtail = (coh*CMh0*Sh/S + (coh*CMh1-xhbox-xCG)*    CLh*Sh/S )/cma
      #     CMfuse = CMVf1*(CL-CLMf0)/(S*cma)

      #---- store results for returning
      #      para[iaCMwing] = CMwing
      #      para[iaCMtail] = CMtail
      #      para[iaCMfuse] = CMfuse

      para[iaxCG] = xCG
      para[iaxCP] = xCP
      para[iaxNP] = xNP

      return
end # balance_aircraft!

"""
    size_htail(ac, paraF, paraB, paraC)

Solves for the feasible horizontal tail area (`Sh`) and wing box location 
(`xwbox`) to ensure: (1) pitch trim req't with forward CG and
(2) stability requirement with aft CG across different flight conditions.

This routine iteratively adjusts:
- Horizontal tail area (`Sh`): Ensuring sufficient control authority.
- Wing box location (`xwbox`): Maintaining static and dynamic stability.

The routine considers:
- Max and min CG locations** (`xcgF`, `xcgB`) computed using `CG_limits()`.
- Aerodynamic parameters (`paraF`, `paraB`, `paraC`).
- Static margin constraints (`SM`).
- Fuel and payload distribution effects.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac` : Aircraft object
    - `paraF[.]` : Aerodynamic parameters for the forward CG trim condition.
    - `paraB[.]` : Aerodynamic parameters for the aft CG stability condition.
    - `paraC[.]` : Aerodynamic parameters for the cruise tail lift condition.

    **Outputs:**  
    No direct return values. Instead, the function updates key fields inside `ac`:
      - `htail.layout.S` : Horizontal tail area (`Sh`).
      - `wing.layout.box_x` : Wing box location (`xwbox`).
      - `wing.layout.x` : Wing location (adjusted for `xwbox`).

**Note**: two flags determine the sizing strategy:

      - `htail.opt_sizing`
            * = "fixed_Vh"    set Sh from prescribed tail volume
            * = "CLmax_fwdCG" set Sh to meet the "worst-case" scenario: max wing CL with most fwd CG, with an assumed max CLh

      - `ac.opt_move_wing`
            * = "fixed"         no changes to wing location
            * = "fixed_CLh"     adjust wingbox location to set tail lift at cruise, CLh = CLhspec
            * = "min_static_margin" adjust wingbox location to get SM = SMmin at aft-CG

The two flags can be set independently and affect how the two stability residuals are driven to zero. The 2x2 system is built sequentially as annotated in the source code for this function. 

"""
function size_htail(ac, paraF, paraB, paraC; Ldebug::Bool = false)
      #TODO find a way to remove the para inputs and use ac instead
      parg, options, fuse, fuse_tank, wing, htail, vtail, engine, landing_gear = unpack_ac_components(ac)

      itmax = 10
      toler = 1.0e-7

      r = zeros(Float64, 2)
      a = zeros(Float64, (2, 2))

      sweep = wing.layout.sweep
      cosL = cosd(sweep)

      #---- set CG limits with worst-case payload arrangements
      rfuelF, rfuelB, rpayF, rpayB, xcgF, xcgB = CG_limits(ac; Ldebug = Ldebug)

      rpayC = 1.0

      parg[igxCGfwd] = xcgF
      parg[igxCGaft] = xcgB

      CLF = paraF[iaclpmax] * cosL^2
      CLhF = htail.CL_max_fwd_CG

      CLC = paraC[iaCL]
      CLhC = parg[igCLhspec]

      SM = htail.SM_min
      cma = wing.mean_aero_chord

      # Unpack Weights
      Wpay = parg[igWpay]
      Wfuel = parg[igWfuel]
      Wfuse = fuse.weight
      Wwing = wing.weight
      Wstrut = wing.strut.weight
      Whtail = htail.weight
      Wvtail = vtail.weight
      Weng = parg[igWeng]

      xWfuel = parg[igxWfuel]

      # Wtshaft = parg[igWtshaft] 
      # Wgen    = parg[igWgen   ] 
      # Winv    = parg[igWinv   ] 
      # Wmot    = parg[igWmot   ] 
      # Wfan    = parg[igWfan   ] 

      Wtesys = parg[igWtesys]
      xWtesys = parg[igxWtesys]

      nftanks = fuse_tank.tank_count
      lftank = parg[iglftank]
      Wftank = parg[igWftank]
      xWftank = parg[igxWftank]

      Whpesys = parg[igWMTO] * fuse.HPE_sys.W
      Wlgnose = landing_gear.nose_gear.weight.W
      Wlgmain = landing_gear.main_gear.weight.W

      xWfuse = fuse.moment

      dxWwing = wing.dxW
      dxWstrut = wing.strut.dxW
      dxWhtail = htail.dxW
      dxWvtail = vtail.dxW

      xeng = parg[igxeng]
      xlgnose = landing_gear.nose_gear.weight.r[1]

      # xtshaft = parg[igxtshaft ]
      # xgen    = parg[igxgen    ]
      # xinv    = parg[igxinv    ]
      # xmot    = parg[igxmot    ]
      # xfan    = parg[igxfan    ]
      xftank = parg[igxftank]

      # Calculate x location of cabin centroid and length of cabin
      xcabin,lcabin = cabin_centroid(nftanks,fuse,parg[igxftankaft],lftank)

      #---- payload CG locations for forward and aft CG locations
      xpayF = xcabin + (0.0 - 0.5) * lcabin * (1.0 - rpayF)
      xpayB = xcabin + (1.0 - 0.5) * lcabin * (1.0 - rpayB)
      xpayC = xcabin

      xwbox = wing.layout.box_x

      #---- wing centroid offset from wingbox, assumed fixed in CG calculations
      dxwing = wing.layout.x - wing.layout.box_x

      #---- main LG offset from wingbox, assumed fixed in CG calculations
      dxlg = parg[igxCGaft] + landing_gear.main_gear.distance_CG_to_landing_gear - wing.layout.box_x

      S = wing.layout.S
      Sh = htail.layout.S
      co = wing.layout.root_chord
      coh = htail.layout.root_chord
      xhbox = htail.layout.box_x
      xvbox = vtail.layout.box_x

      WfuelC = paraC[iafracW] * parg[igWMTO] -
               rpayC * Wpay -
               Wfuse -
               Wtesys -
               Wftank -
               Wwing -
               Wstrut -
               Whtail -
               Wvtail -
               Weng -
               Whpesys -
               Wlgnose -
               Wlgmain

      rfuelC = WfuelC / Wfuel


      #---- for better convergence, 
      #      will assume that HT weight and moment will scale with its area
      Sh_o = Sh
      Whtail_o = Whtail
      dxWhtail_o = dxWhtail

      #c      do 100 iter = 0, -2, -1

      xWF, xWB, xWc = 0.0, 0.0, 0.0
      WF, WB, Wc = 0.0, 0.0, 0.0
      dmax = 0.0

      @inbounds for iter = 1:itmax

            # HT weight, moment and derivatives ∂Whtail/∂Sh
            Whtail = (Whtail_o / Sh_o) * Sh
            Whtail_Sh = Whtail_o / Sh_o

            dxWhtail = (dxWhtail_o / Sh_o) * Sh
            dxWhtail_Sh = dxWhtail_o / Sh_o

            #---- empty (no fuel, no payload) weight
            We = Wfuse +
                 Wwing +
                 Wtesys +
                 Wftank +
                 Wstrut +
                 Whtail +
                 Wvtail +
                 Weng +
                 Whpesys +
                 Wlgnose +
                 Wlgmain

            We_Sh = Whtail_Sh

            #---- empty (no-payload) weight moment
            xWe = xWfuse + xWtesys + xWftank +
                  Wwing * xwbox + dxWwing +
                  Wstrut * xwbox + dxWstrut +
                  Whtail * xhbox + dxWhtail +
                  Wvtail * xvbox + dxWvtail +
                  Weng * xeng +
                  Whpesys * fuse.HPE_sys.r.x +
                  landing_gear.nose_gear.moment +
                  Wlgmain * (xwbox + dxlg)

            xWe_Sh = Whtail_Sh * xhbox + dxWhtail_Sh
            xWe_xw = Wwing +
                     Wstrut +
                     Wlgmain


            #---- total weight at forward and aft CG limits, and cruise
            WF = We + rfuelF * Wfuel + rpayF * Wpay
            WB = We + rfuelB * Wfuel + rpayB * Wpay
            WC = We + rfuelC * Wfuel + rpayC * Wpay

            WF_Sh = We_Sh
            WB_Sh = We_Sh
            WC_Sh = We_Sh

            #---- total weight moment at forward and aft CG limits
            if !(ac.options.has_wing_fuel)
                  xWfuel_xw = 0.0
            else
                  xWfuel_xw = Wfuel
            end
            xWF = xWe + rpayF * Wpay * xpayF + rfuelF * xWfuel
            xWB = xWe + rpayB * Wpay * xpayB + rfuelB * xWfuel
            xWC = xWe + rpayC * Wpay * xpayC + rfuelC * xWfuel

            xWF_xw = xWe_xw + rfuelF * xWfuel_xw
            xWB_xw = xWe_xw + rfuelB * xWfuel_xw
            xWC_xw = xWe_xw + rfuelC * xWfuel_xw

            xWF_Sh = xWe_Sh
            xWB_Sh = xWe_Sh
            xWC_Sh = xWe_Sh

            if compare_strings(htail.opt_sizing, "fixed_Vh")
                  #----- set HT area from fixed volume (Section 2.12.1 of TASOPT docs)
                  lhtail = htail.layout.x - (xwbox + dxwing)
                  lhtail_xw = -1.0

                  Vh = htail.volume     # Tail volume is specified here!
                  Sh = Vh * S * cma / lhtail # Corresponding tail surf area

                  r[1] = Sh * lhtail - Vh * S * cma  # This is 0 (residual is 0 ∵ Vh = Sh*lh/(S*cma))
                  a[1, 1] = lhtail
                  a[1, 2] = Sh * lhtail_xw

            elseif compare_strings(htail.opt_sizing, "CLmax_fwdCG")
                  #----- set HT area by worst case: max pitch trim power at forward CG case w/ max wing CL
                  CMw0 = paraF[iaCMw0]
                  CMw1 = paraF[iaCMw1]
                  CMh0 = paraF[iaCMh0]
                  CMh1 = paraF[iaCMh1]
                  CMVf1 = parg[igCMVf1]
                  CLMf0 = parg[igCLMf0]
                  CL = CLF
                  CLh = CLhF

                  cCM = co * CMw0 + (co * CMw1 - xwbox) * (CL - CLh * Sh / S) +
                        coh * CMh0 * Sh / S + (coh * CMh1 - xhbox) * CLh * Sh / S +
                        CMVf1 * (CL - CLMf0) / S
                  cCM_xw = -(CL - CLh * Sh / S)
                  cCM_Sh = (co * CMw1 - xwbox) * (-CLh / S) +
                           coh * CMh0 / S + (coh * CMh1 - xhbox) * CLh / S

                  # Residual eqn ℛₘ that is equivalent to CP and CG coinciding
                  r[1] = cCM / CL + xWF / WF
                  a[1, 1] = cCM_Sh / CL + xWF_Sh / WF - (xWF / WF^2) * WF_Sh  # ∂ℛₘ/∂Sₕ
                  a[1, 2] = cCM_xw / CL + xWF_xw / WF                     # ∂ℛₘ/∂xw

                  # +--                 --+  +-   -+      +-   -+
                  # | ∂ℛₘ/∂Sₕ    ∂ℛₘ/∂xw |  | δSₕ  |      | ℛₘ |
                  # |                     |  |     |  =   |     |
                  # |                     |  | δxw |      |     |
                  # +--                 --+  +-   -+      +-   -+

            end


            if compare_strings(options.opt_move_wing, "fixed")
                  #----- fix wing location
                  r[2] = 0.0
                  a[2, 1] = 0.0
                  a[2, 2] = 1.0

            elseif compare_strings(options.opt_move_wing, "fixed_CLh")
                  #----- set wing location to get CLh=CLhspec in cruise
                  CMw0 = paraC[iaCMw0]
                  CMw1 = paraC[iaCMw1]
                  CMh0 = paraC[iaCMh0]
                  CMh1 = paraC[iaCMh1]
                  CMVf1 = parg[igCMVf1]
                  CLMf0 = parg[igCLMf0]
                  CL = CLC
                  CLh = CLhC  #CLh is specified here

                  cCM = co * CMw0 + (co * CMw1 - xwbox) * (CL - CLh * Sh / S) +
                        coh * CMh0 * Sh / S + (coh * CMh1 - xhbox) * CLh * Sh / S +
                        CMVf1 * (CL - CLMf0) / S
                  #c   &        - neng*xenginl*CLn*Afan/S
                  cCM_xw = -(CL - CLh * Sh / S)
                  cCM_Sh = (co * CMw1 - xwbox) * (-CLh / S) +
                           coh * CMh0 / S + (coh * CMh1 - xhbox) * CLh / S

                  r[2] = cCM / CL + xWC / WC
                  a[2, 1] = cCM_Sh / CL + xWC_Sh / WC - (xWC / WC^2) * WC_Sh
                  a[2, 2] = cCM_xw / CL + xWC_xw / WC

            elseif compare_strings(options.opt_move_wing, "min_static_margin")
                  #----- set wing location by stability margin at aft-CG case
                  CMw1 = paraB[iaCMw1]
                  CMh1 = paraB[iaCMh1]
                  CMVf1 = parg[igCMVf1]
                  dCLhdCL = parg[igdCLhdCL]

                  dCLndCL = parg[igdCLndCL]
                  neng = parg[igneng]
                  xeng = parg[igxeng]
                  dfan = parg[igdfan]
                  Afan = 0.25 * pi * dfan^2
                  xengcp = xeng - 0.25 * dfan
                  # Eqn 297, 298:
                  xNP = (co * CMw1 - xwbox) * (dCLhdCL * Sh / S - 1.0) -
                        (coh * CMh1 - xhbox) * dCLhdCL * Sh / S -
                        CMVf1 / S -
                        neng * xengcp * dCLndCL * Afan / S  #∂CLnace/∂CL moments caused by engines

                  xNP_xw = -(dCLhdCL * Sh / S - 1.0)
                  xNP_Sh = (co * CMw1 - xwbox) * dCLhdCL / S -
                           (coh * CMh1 - xhbox) * dCLhdCL / S

                  # ℛsh(xw, Sh) ≡ xCG - xNP + fSM*cma 
                  # → must have some min static pitch stability i.e. xCG in front of xNP by a static margin
                  r[2] = xWB / WB - xNP + SM * cma
                  a[2, 1] = xWB_Sh / WB - xNP_Sh - (xWB / WB^2) * WB_Sh
                  a[2, 2] = xWB_xw / WB - xNP_xw

                  # +-                   -+  +-   -+      +-   -+
                  # | ∂ℛₘ/∂Sₕ    ∂ℛₘ/∂xw |  | δSₕ  |      | ℛₘ |
                  # |                     |  |     |  =   |     |
                  # | ∂ℛsh/∂Sₕ   ∂ℛsh/∂xw|  | δxw |      |ℛsh |
                  # +-                   -+  +-   -+      +-   -+

            end

            det = a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1]
            d1 = (r[1] * a[2, 2] - a[1, 2] * r[2]) / det
            d2 = (a[1, 1] * r[2] - r[1] * a[2, 1]) / det

            dSh = -d1
            dxw = -d2

            Sh = Sh + dSh
            xwbox = xwbox + dxw

            dmax = max(abs(dxw) / cma, abs(dSh) / (0.2 * S))
            if (dmax < toler)
                  break
            end

      end #end for loop iter
      if (dmax > toler)
            @warn "`size_htail()`: Pitch not converged. dxwbox,dSh = $dxw, $dSh"
      end


      #---- set converged results
      htail.layout.S = Sh
      wing.layout.box_x = xwbox
      wing.layout.x = xwbox + dxwing

      #Move engine as well to maintain the input offset distance from engine to wing box
      dxeng2wbox = parg[igdxeng2wbox]
      parg[igxeng] = xwbox - dxeng2wbox

      if compare_strings(ac.htail.opt_sizing, "fixed_Vh")
            #----- for fixed HT area, find minimum required CLh for foward-CG trim
            CMw0 = paraF[iaCMw0]
            CMw1 = paraF[iaCMw1]
            CMh0 = paraF[iaCMh0]
            CMh1 = paraF[iaCMh1]
            CMVf1 = parg[igCMVf1]
            CLMf0 = parg[igCLMf0]
            CL = CLF
            CLh = CLhF
            cCM = co * CMw0 + (co * CMw1 - xwbox) * (CL - CLh * Sh / S) +
                  coh * CMh0 * Sh / S + (coh * CMh1 - xhbox) * CLh * Sh / S +
                  CMVf1 * (CL - CLMf0) / S
            cCM_CLh = (co * CMw1 - xwbox) * (-Sh / S) +
                      (coh * CMh1 - xhbox) * Sh / S

            res = cCM / CL + xWF / WF
            res_CLh = cCM_CLh / CL
            dCLh = -res / res_CLh
            CLh = CLh + dCLh
            htail.CL_max_fwd_CG = CLh
      end


      return
end # size_htail

"""
    CG_limits(ac; Ldebug::Bool = false)

Computes the most forward (`xcgF`) and most rearward (`xcgB`) 
center of gravity (CG) locations based on payload extremes,
along with the corresponding payload fractions. Formerly, `cglpay()`.

## Description
This function determines the CG shift due to varying passenger and fuel load configurations.
- `xcgF`: The forward-most CG position, assuming worst-case forward payload arrangement.
- `xcgB`: The rearward-most CG position, assuming worst-case aft payload arrangement.
- `rpayF`: The fraction of payload contributing to the forward-most CG.
- `rpayB`: The fraction of payload contributing to the rearward-most CG.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac` : Aircraft object
    - `Ldebug` : Debug flag for verbose output (default: false, optional)

    **Outputs:**
    Returns six values:
    - `rfuelF` : Fuel fraction for forward CG case (always 0.0).
    - `rfuelB` : Fuel fraction for aft CG case (always 0.0).
    - `rpayF` : Payload fraction contributing to the forward CG (`xcgF`).
    - `rpayB` : Payload fraction contributing to the rearward CG (`xcgB`).
    - `xcgF` : Most forward CG location.
    - `xcgB` : Most rearward CG location.

"""
function CG_limits(ac; Ldebug::Bool = false)
      parg, options, fuse, fuse_tank, wing, htail, vtail, engine, landing_gear = unpack_ac_components(ac)

      Wpay = parg[igWpay]
      Wfuel = parg[igWfuel]
      Wfuse = fuse.weight
      Wwing = wing.weight
      Wstrut = wing.strut.weight
      Whtail = htail.weight
      Wvtail = vtail.weight
      Weng = parg[igWeng]

      xWfuel = parg[igxWfuel]

      Wtesys = parg[igWtesys]
      #      xWtesys = parg[igxWtesys]

      nftanks = fuse_tank.tank_count
      lftank = parg[iglftank]
      Wftank = parg[igWftank]
      #      xWftank = parg[igxWftank]

      Whpesys = parg[igWMTO] * fuse.HPE_sys.W
      Wlgnose = landing_gear.nose_gear.weight.W
      Wlgmain = landing_gear.main_gear.weight.W

      xcabin,lcabin = cabin_centroid(nftanks,fuse,parg[igxftankaft],lftank)
      delxw = wing.layout.x - wing.layout.box_x

      #---- zero fuel is assumed to be worst case forward and full fuel worst case aft
      rfuel = 0.0

      # See Eqn 283 in TASOPT documentation

      We = rfuel * Wfuel +
           Wfuse +
           Wtesys +
           Wftank +
           Wwing +
           Wstrut +
           Whtail +
           Wvtail +
           Weng +
           Whpesys +
           Wlgnose +
           Wlgmain

      xWe = rfuel * xWfuel +
            fuse.moment + parg[igxWtesys] + parg[igxWftank] +
            Wwing * wing.layout.box_x + wing.dxW +
            Wstrut * wing.layout.box_x + wing.strut.dxW +
            Whtail * htail.layout.box_x + htail.dxW +
            Wvtail * vtail.layout.box_x + vtail.dxW +
            Weng * parg[igxeng] +
            Whpesys * fuse.HPE_sys.r.x +
            landing_gear.nose_gear.moment +
            Wlgmain * (wing.layout.box_x + delxw + landing_gear.main_gear.distance_CG_to_landing_gear)


      # Some derivation here:        
      #   xpay = xcabin + (ξpay-0.5)*lcabin*(1.0-rpay)  [1]                                   
      #    W =   We +      Wpay*rpay + rf* Wf           [2]                         
      #   xW =  xWe + xpay*Wpay*rpay + rf*xWf           [3]                         
      #  xcg = xW/W                                     [4]
      #
      #  Now, xpay*Wpay*rpay = (xcabin + (ξpay-0.5)*lcabin)*Wpay*rpay - (ξpay-0.5)*lcabin*Wpay*rpay²
      #                      =               a₁                 *rpay +            a₂         *rpay²
      #
      #  ∴ xW = a₀ + a₁rpay + a₂rpay² ;  a₀ = xWe + rf*xWf
      #     W = b₀ + b₁rpay           ;  b₀ =  We + rf* Wf, b₁ = Wpay

      # Explanation:
      #  now we find the stationary points of 
      #  xcg wrt to rpay. xcg = xW/W by setting 
      #  ∂xcg/∂rpay = 0
      #  this simplifies to a quadratic equation that gives rpay for the 
      #  min and max values of xcg
      #
      #  Repeat this for ξ = 0 and ξ = 1 (passengers all in front vs all in the back)
      #
      ξ = [0.0, 1.0]
      sgn = [-1.0, 1.0]

      xftank = parg[igxWftank] / Wftank
      if (options.has_wing_fuel) #Fuel is in wings
            rf = [0.0, 0.0]
      elseif xftank < xcabin
            rf = [1.0, 0.0]
      else
            rf = [0.0, 1.0]
      end
      
      rpay = zeros(Float64, 2)
      xcg = zeros(Float64, 2)

      @inbounds for i = 1:2

            a0 = xWe + rf[i] * xWfuel
            b0 = We + rf[i] * Wfuel
            b1 = Wpay

            a1 = (xcabin + (ξ[i] - 0.5) * lcabin) * Wpay
            a2 = -(ξ[i] - 0.5) * lcabin * Wpay

            AA = a2 * b1
            BB = 2.0 * a2 * b0
            CC = a1 * b0 - a0 * b1

            if BB^2 - 4.0 * AA * CC ≤ 0.0 && Ldebug
                  @warn "CG_limits(): Warning: No real roots for quadratic equation"
                  println("a2 = $a2 ; b0^2 = $(b0^2); a0 = $(a0); b1^2 = $(b1^2); a1 = $(a1)")
            end

            rpay[i] = (-BB - sgn[i] * sqrt(BB^2 - 4.0 * AA * CC)) * 0.5 / AA

            xpay = xcabin + (ξ[i] - 0.5) * lcabin * (1.0 - rpay[i])
            W = Wpay * rpay[i] + We + rf[i] * Wfuel
            xW = xpay * Wpay * rpay[i] + xWe + rf[i] * xWfuel

            xcg[i] = xW / W
      end

      #      rpayF = rpay[1]
      #      rpayB = rpay[2]
      #
      #      xcgF = xcg[1]
      #      xcgB = xcg[2]

      return rf[1], rf[2], rpay[1], rpay[2], xcg[1], xcg[2]
end # CG_limits



"""
    cabin_centroid(nftanks, fuse, xftankaft, lftank)

Computes the centroid (`xcabin`) and length (`lcabin`) of the passenger cabin
accounting for the presence and location of fuel tanks.

determines the cabin centroid (`xcabin`) and cabin length (`lcabin`) based on:
- The number of fuel tanks (`nftanks`).
- Whether the fuel tank is located at the front or rear (`xftankaft`).
- The length of the fuel tank (`lftank`).
- The fuselage layout (`fuse.layout`), including pressure shell and cylindrical section dimensions.

The cabin centroid is calculated as the midpoint of the effective passenger cabin length,
which varies depending on fuel tank placement.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `nftanks` : Number of fuel tanks (`0`, `1`, or `2`).
    - `fuse` : Fuselage object containing geometry properties.
    - `xftankaft` : Binary flag (`0.0` if fuel tank is at the front, otherwise rear).
    - `lftank` : Length of the fuel tank.

    **Outputs:**
    - `xcabin` : `x`-coordinate of the cabin centroid.
    - `lcabin` : Length of the passenger cabin.
"""
function cabin_centroid(nftanks,fuse,xftankaft,lftank)
      # Calculate x location of cabin centroid  and length of cabin
      if nftanks == 1
            if xftankaft == 0.0 #If tank is at the front
                  xcabin = 0.5 * (fuse.layout.x_start_cylinder + lftank + 2.0*ft_to_m + fuse.layout.x_pressure_shell_aft)
                  lcabin = fuse.layout.x_pressure_shell_aft - fuse.layout.x_end_cylinder + fuse.layout.l_cabin_cylinder #cabin length is smaller if there are fuel tanks
            else #tank is at rear
                  xcabin = 0.5 * (fuse.layout.x_pressure_shell_fwd + fuse.layout.x_end_cylinder - (lftank + 2.0*ft_to_m))
                  lcabin = fuse.layout.x_start_cylinder - fuse.layout.x_pressure_shell_fwd + fuse.layout.l_cabin_cylinder #cabin length is smaller if there are fuel tanks
            end
      elseif nftanks == 2
            xcabin = 0.5 * (fuse.layout.x_pressure_shell_fwd + fuse.layout.x_pressure_shell_aft) #TODO noticed convergence issues if the average of the blends is used instead
            lcabin = fuse.layout.l_cabin_cylinder
      elseif nftanks == 0
            xcabin = 0.5 * (fuse.layout.x_pressure_shell_fwd + fuse.layout.x_pressure_shell_aft)
            lcabin = fuse.layout.x_pressure_shell_aft - fuse.layout.x_pressure_shell_fwd
      end
      return xcabin,lcabin
end