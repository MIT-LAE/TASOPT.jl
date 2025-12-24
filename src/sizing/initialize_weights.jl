"""
    initialize_sizing_loop!(ac)

Initializes weights, geometry estimates, and parameters for the first iteration
of the aircraft sizing loop. Sets all state directly in the aircraft structs and
parameter arrays.

This function computes initial estimates using:
- Payload-based weight fractions for components
- Breguet range equation for fuel fraction
- Simple geometric relationships for wing/tail sizing

All values are approximate and will be refined during iteration.
"""
function initialize_sizing_loop!(ac)
    # Unpack aircraft components and parameter arrays
    imission = 1
    parg, parm, para, pare, _, fuse, _, wing, htail, vtail, _, landing_gear = unpack_ac(ac, imission)

    # Extract design parameters
    Wpay = parm[imWpay]
    Wpaymax = parg[igWpaymax]
    Rangetot = parm[imRange]
    freserve = parg[igfreserve]
    neng = parg[igneng]

    # Extract layout parameters
    xhbox = htail.layout.box_x
    xvbox = vtail.layout.box_x

    # Landing gear weight fractions
    flgnose = landing_gear.nose_gear.overall_mass_fraction
    flgmain = landing_gear.main_gear.overall_mass_fraction

    # Dynamic pressure at never-exceed speed
    Vne = parg[igVne]
    qne = 0.5 * ρSL * Vne^2

    # Structural weight factor
    sigfac = parg[igsigfac]

    # ===== Initial weight estimates =====
    Whtail = 0.05 * Wpay / sigfac
    Wvtail = Whtail
    Wwing = 0.5 * Wpay / sigfac
    Wstrut = 0.0
    Wftank = 0.0
    Weng = 0.0

    # ===== Wing geometry estimates =====
    ip = ipcruise1
    W_estimate = 5.0 * Wpay  # Initial weight estimate
    S = W_estimate / (0.5 * pare[ierho0, ip] * pare[ieu0, ip]^2 * para[iaCL, ip])
    b = sqrt(S * wing.layout.AR)
    bs = b * wing.layout.ηs

    # Wing panel weights
    Winn = 0.15 * Wpay / sigfac
    Wout = 0.05 * Wpay / sigfac
    dyWinn = Winn * 0.30 * (0.5 * (bs - wing.layout.root_span))
    dyWout = Wout * 0.25 * (0.5 * (b - bs))

    # ===== Store initial weights in components =====
    htail.weight = Whtail
    vtail.weight = Wvtail
    wing.weight = Wwing
    wing.strut.weight = Wstrut
    wing.inboard.weight = Winn
    wing.outboard.weight = Wout

    parg[igWeng] = Weng
    parg[igWftank] = Wftank

    # Weight moments
    htail.dxW = 0.0
    vtail.dxW = 0.0
    wing.inboard.dyW = dyWinn
    wing.outboard.dyW = dyWout

    # Wing centroid x-offset from wingbox
    calculate_centroid_offset!(wing, b=b, bs=bs)

    # Tail area centroid locations
    htail.layout.x = xhbox
    vtail.layout.x = xvbox

    # ===== Nacelle parameters =====
    parg[iglnace] = 0.5 * S / b
    parg[igfSnace] = 0.2

    # ===== Fuel fraction from Breguet Range Equation =====
    LoD = 18.0  # Initial L/D estimate
    TSFC = 1.0 / 7000.0  # Typical TSFC
    V = pare[ieu0, ipcruise1]
    ffburn = min((1.0 - exp(-Rangetot * TSFC / (V * LoD))), 0.8 / (1.0 + freserve))

    # Mission-point fuel fractions
    ffuelb = ffburn * (1.0 + freserve)   # Start of climb
    ffuelc = ffburn * (0.90 + freserve)  # Start of cruise
    ffueld = ffburn * (0.02 + freserve)  # Start of descent
    ffuele = ffburn * (0.0 + freserve)   # End of descent (landing)

    ffuel = ffuelb  # Max fuel fraction at start of climb

    # Set initial climb γ = 0
    para[iagamV, :] .= 0.0

    # Set initial weight fractions for takeoff/cutback
    para[iafracW, ipstatic:ipcutback] .= 1.0

    # Interpolate weight fractions for climb, cruise, and descent
    interp_Wfrac!(para, ipclimb1, ipclimbn, ffuelb, ffuelc, iafracW, ffuel)
    interp_Wfrac!(para, ipcruise1, ipcruisen, ffuelc, ffueld, iafracW, ffuel)
    interp_Wfrac!(para, ipdescent1, ipdescentn, ffueld, ffuele, iafracW, ffuel)

    # ===== Initial tail areas =====
    Sh = (2.0 * Wpaymax) / (qne * htail.CL_max)
    Sv = (2.0 * Wpaymax) / (qne * vtail.CL_max)

    htail.layout.S = Sh
    vtail.layout.S = Sv

    # ===== Initialize pitch moments =====
    para[iaCMw0:iaCMh1, :] .= 0.0
    para[iaCLh, :] .= 0.0

    # ===== Cruise-climb parameters =====
    gamVcr = 0.0002
    para[iaCD, ipcruise1] = para[iaCL, ipcruise1] / LoD
    para[iagamV, ipcruise1] = gamVcr

    # End-of-cruise pressure
    p0c = pare[iep0, ipcruise1]
    p0d = p0c * (1.0 - ffuel + ffueld) / (1.0 - ffuel + ffuelc)
    pare[iep0, ipcruisen] = p0d

    # ===== Initial OEI thrust and fan sizing =====
    pare[ieFe, iprotate] = 2.0 * Wpay
    pare[ieu0, iprotate] = 70.0
    Afan = 3.0e-5 * Wpay / neng
    parg[igdfan] = sqrt(Afan * 4.0 / π)

    # Fan face Mach numbers
    M2des = 0.6
    pare[ieM2, ipstatic:ipcruisen] .= M2des
    pare[ieM2, ipdescent1:ipdescentn] .= 0.8 * M2des

    # ===== Cooling mass flow ratios =====
    ip = iprotate
    cpc, cp4 = 1080.0, 1340.0
    Rgc, Rg4 = 288.0, 288.0

    M0to = pare[ieu0, ip] / pare[iea0, ip]
    T0to = pare[ieT0, ip]
    epolhc = pare[ieepolhc, ip]
    OPRto = pare[iepilc, ipcruise1] * pare[iepihc, ipcruise1]
    Tt4to = pare[ieTt4, ip]
    dTstrk = pare[iedTstrk, ip]
    Mtexit = pare[ieMtexit, ip]
    efilm = pare[ieefilm, ip]
    tfilm = pare[ietfilm, ip]
    StA = pare[ieStA, ip]

    # Create working arrays for cooling calculation
    Tmrow = fill(parg[igTmetal], ncrowx)

    Tt2to = T0to * (1.0 + 0.5 * (gamSL - 1.0) * M0to^2)
    Tt3to = Tt2to * OPRto^(Rgc / (epolhc * cpc))
    Trrat = 1.0 / (1.0 + 0.5 * Rg4 / (cp4 - Rg4) * Mtexit^2)

    ncrow, epsrow, _, _, _ = mcool(ncrowx, Tmrow,
        Tt3to, Tt4to, dTstrk, Trrat, efilm, tfilm, StA)

    epstot = sum(epsrow[1:ncrow])
    fo = pare[iemofft, ip] / pare[iemcore, ip]
    fc = (1.0 - fo) * epstot

    # Store cooling parameters for all mission points
    for jp in 1:iptotal
        pare[iefc, jp] = fc
        pare[ieepsc1:ieepsc1+ncrowx-1, jp] .= epsrow
        pare[ieTmet1:ieTmet1+ncrowx-1, jp] .= Tmrow
    end

    # ===== Initial WMTO estimate =====
    # Use assumed engine fraction and Breguet fuel fraction
    feng = 0.08
    fsum = feng + ffuel + fuse.HPE_sys.W + flgnose + flgmain

    # Compute initial WMTO from fixed weights and fractions
    WMTO = (Wpay + fuse.weight + Wwing + Wstrut + Whtail + Wvtail) / (1.0 - fsum)

    # Store initial estimates
    parg[igWMTO] = WMTO
    parg[igWeng] = WMTO * feng
    parg[igWfuel] = WMTO * ffuel

    return nothing
end

"""
    interp_Wfrac!(para, ip_start, ip_end, ffuel1, ffuel2, iafracW, ffuel)

Interpolates iafracW from two mission points
"""
function interp_Wfrac!(para, ip_start, ip_end, ffuel1, ffuel2, iafracW, ffuel)
    @inbounds for ip in ip_start:ip_end
        frac = float(ip - ip_start) / float(ip_end - ip_start)
        ffp = ffuel1 * (1.0 - frac) + ffuel2 * frac
        para[iafracW, ip] = 1.0 - ffuel + ffp
    end
end