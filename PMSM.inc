# ------------- MOTOR ----------------
# Parameters for the Turbo-electric propulsion system
parmot = zeros(Float64, imot_total)
# PMSM parameters
parmot[imot_ratSM  ]  = 0.86    # Slot to magnet height ratio (hS/hM)
parmot[imot_ratAg  ]  = 1/160   # Air-gap ratio (hAg - hRS)/dRot
parmot[imot_ratM   ]  = 6       # Magnet height ratio hM/hAg
parmot[imot_ratSp  ]  = 36/6    # Number of slots to pole-pairs NS/p
parmot[imot_ratW   ]  = 0.5     # Ratio of tooth width to slot pitch = wT/wST
parmot[imot_ratShft]  = 3/4     # Ratio of inner to outer dia for hollow Shaft
parmot[imot_z      ]  = 3       # Phases

# parmot[imot_kServ  ]  = 1.15    # Rucker2005
# parmot[imot_kServ  ]  = 1.5     # Yoon2016 (Tab.4) Excluding Ground Cylinder and Heat Sink
parmot[imot_kServ  ]  = 1.7     # Yoon2016 (Tab.4) Including Ground Cylinder and Heat Sink

parmot[imot_ratSd ]  = 1/50    # Slot depression to height ratio hsd/hsd
parmot[imot_kpf   ]  = 0.5     # kpf = syrat = Packing factor
parmot[imot_wSd   ]  = 1e-6    # Slot depression width
parmot[imot_Nshrt ]  = 1       # Number of short pitched slots
parmot[imot_hRS   ]  = 0.002   # [m] Thickness of retaining sleeve for Vtip = 200 m/s
parmot[imot_Tarm  ]  = 273.16 + 80 # Armature Temperature in K
parmot[imot_kst   ]  = 0.98    # Stacking factor for silicon electrical steel

parmot[imot_p]     = 8   # pole-pairs
        
parmot[imot_Br  ]  = 1.2 # Magnet remanence (Nd-Fe-B magnets)
parmot[imot_BSat]  = 1.6 # Can go upto 1.9 perhaps
parmot[imot_mu0 ]  = 4π * 1e-7 # Vacuum permeability

parmot[imot_rhoMag  ]     = 7400. # Ne-Fe-B magnet density [kg/m^3]
parmot[imot_rhoCu   ]     = 8900. # Cu density
parmot[imot_rhoFe   ]     = 7700. # Silicon electric Iron density
parmot[imot_rhoSteel]     = 7850. # S275 steel 
parmot[imot_tauMax  ]     = 110.0e6 # Max shaft torsion shear stress S275 Steel
parmot[imot_sigCu   ]     = 1.0e7   # Electrical condutivity of copper at 293K
parmot[imot_thetaCu ]     = 3.81e-3 # Temp coeff of copper (Hanselman)

parmot[imot_psi]     = 0.0 # Power factor angle [rad]

parmot[imot_pb0 ]   = 36.79  # Base factor for power
parmot[imot_Bb0 ]   = 1.0    # Base factor for flux density
parmot[imot_fb0 ]   = 1000.  # Base factor for frequency
parmot[imot_epsb]   = 2.12   # Base factor for flux density exponent
parmot[imot_epsf]   = 1.68   # Base factor for frequency exponent


# ------------- GENERATOR ----------------
# Parameters for the Turbo-electric propulsion system
pargen = zeros(Float64, igen_total)
# PMSM parameters
pargen[igen_ratSM  ]  = 0.86    # Slot to magnet height ratio (hS/hM)
pargen[igen_ratAg  ]  = 1/160   # Air-gap ratio (hAg - hRS)/dRot
pargen[igen_ratM   ]  = 6       # Magnet height ratio hM/hAg
pargen[igen_ratSp  ]  = 36/6    # Number of slots to pole-pairs NS/p
pargen[igen_ratW   ]  = 0.5     # Ratio of tooth width to slot pitch = wT/wST
pargen[igen_ratShft]  = 3/4     # Ratio of inner to outer dia for hollow Shaft
pargen[igen_z      ]  = 3       # Phases

# pargen[igen_kServ  ]  = 1.15    # Rucker2005
# pargen[igen_kServ  ]  = 1.5     # Yoon2016 (Tab.4) Excluding Ground Cylinder and Heat Sink
pargen[igen_kServ  ]  = 1.7     # Yoon2016 (Tab.4) Including Ground Cylinder and Heat Sink

pargen[igen_ratSd ]  = 1/50    # Slot depression to height ratio hsd/hsd
pargen[igen_kpf   ]  = 0.5     # kpf = syrat = Packing factor
pargen[igen_wSd   ]  = 1e-6    # Slot depression width
pargen[igen_Nshrt ]  = 1       # Number of short pitched slots
pargen[igen_hRS   ]  = 0.002   # [m] Thickness of retaining sleeve for Vtip = 200 m/s
pargen[igen_Tarm  ]  = 273.16 + 80 # Armature Temperature in K
pargen[igen_kst   ]  = 0.98    # Stacking factor for silicon electrical steel

pargen[igen_p]     = 8   # pole-pairs
        
pargen[igen_Br  ]  = 1.2 # Magnet remanence (Nd-Fe-B magnets)
pargen[igen_BSat]  = 1.6 # Can go upto 1.9 perhaps
pargen[igen_mu0 ]  = 4π * 1e-7 # Vacuum permeability

pargen[igen_rhoMag  ]     = 7400. # Ne-Fe-B magnet density [kg/m^3]
pargen[igen_rhoCu   ]     = 8900. # Cu density
pargen[igen_rhoFe   ]     = 7700. # Silicon electric Iron density
pargen[igen_rhoSteel]     = 7850. # S275 steel 
pargen[igen_tauMax  ]     = 110.0e6 # Max shaft torsion shear stress S275 Steel
pargen[igen_sigCu   ]     = 1.0e7   # Electrical condutivity of copper at 293K
pargen[igen_thetaCu ]     = 3.81e-3 # Temp coeff of copper (Hanselman)

pargen[igen_psi]     = 0.0 # Power factor angle [rad]

pargen[igen_pb0 ]   = 36.79  # Base factor for power
pargen[igen_Bb0 ]   = 1.0    # Base factor for flux density
pargen[igen_fb0 ]   = 1000.  # Base factor for frequency
pargen[igen_epsb]   = 2.12   # Base factor for flux density exponent
pargen[igen_epsf]   = 1.68   # Base factor for frequency exponent


