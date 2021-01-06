
# Parameters for the Turbo-electric propulsion system
parte = zeros(Float64, ite_total)
# PMSM parameters
parte[ite_ratSM  ]  = 0.86    # Slot to magnet height ratio (hS/hM)
parte[ite_ratAg  ]  = 1/160   # Air-gap ratio (hAg - hRS)/dRot
parte[ite_ratM   ]  = 6       # Magnet height ratio hM/hAg
parte[ite_ratSp  ]  = 36/6    # Number of slots to pole-pairs NS/p
parte[ite_ratW   ]  = 0.5     # Ratio of tooth width to slot pitch = wT/wST
parte[ite_ratShft]  = 3/4     # Ratio of inner to outer dia for hollow Shaft
parte[ite_z      ]  = 3       # Phases

# parte[ite_kServ  ]  = 1.15    # Rucker2005
# parte[ite_kServ  ]  = 1.5     # Yoon2016 (Tab.4) Excluding Ground Cylinder and Heat Sink
parte[ite_kServ  ]  = 1.7     # Yoon2016 (Tab.4) Including Ground Cylinder and Heat Sink

parte[ite_ratSd ]  = 1/50    # Slot depression to height ratio hsd/hsd
parte[ite_kpf   ]  = 0.5     # kpf = syrat = Packing factor
parte[ite_wSd   ]  = 1e-6    # Slot depression width
parte[ite_Nshrt ]  = 1       # Number of short pitched slots
parte[ite_hRS   ]  = 0.002   # [m] Thickness of retaining sleeve for Vtip = 200 m/s
parte[ite_Tarm  ]  = 273.16 + 80 # Armature Temperature in K
parte[ite_kst   ]  = 0.98    # Stacking factor for silicon electrical steel

parte[ite_p]     = 8   # pole-pairs
        
parte[ite_Br  ]  = 1.2 # Magnet remanence (Nd-Fe-B magnets)
parte[ite_BSat]  = 1.6 # Can go upto 1.9 perhaps
parte[ite_mu0 ]  = 4π * 1e-7 # Vacuum permeability

parte[ite_rhoMag  ]     = 7400. # Ne-Fe-B magnet density [kg/m^3]
parte[ite_rhoCu   ]     = 8900. # Cu density
parte[ite_rhoFe   ]     = 7700. # Silicon electric Iron density
parte[ite_rhoSteel]     = 7850. # S275 steel 
parte[ite_tauMax  ]     = 110.0e6 # Max shaft torsion shear stress S275 Steel
parte[ite_sigCu   ]     = 1.0e7   # Electrical condutivity of copper at 293K
parte[ite_thetaCu ]     = 3.81e-3 # Temp coeff of copper (Hanselman)

parte[ite_psi]     = 0.0 # Power factor angle [rad]

parte[ite_pb0 ]   = 36.79  # Base factor for power
parte[ite_Bb0 ]   = 1.0    # Base factor for flux density
parte[ite_fb0 ]   = 1000.  # Base factor for frequency
parte[ite_epsb]   = 2.12   # Base factor for flux density exponent
parte[ite_epsf]   = 1.68   # Base factor for frequency exponent


