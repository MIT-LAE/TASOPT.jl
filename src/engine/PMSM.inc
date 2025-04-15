# WARNING: this file format is deprecated for TASOPT.jl (no turboelectric par array support)
# File is retained in case it's useful for future TE development

# ------------- MOTOR ----------------
# Parameters for the Turbo-electric propulsion system
parmot = zeros(Float64, ite_total)
# PMSM parameters
parmot[ite_ratSM  ]  = 0.86    # Slot to magnet height ratio (hS/hM)
parmot[ite_ratAg  ]  = 1/160   # Air-gap ratio (hAg - hRS)/dRot
parmot[ite_ratM   ]  = 6       # Magnet height ratio hM/hAg
parmot[ite_ratSp  ]  = 36/6    # Number of slots to pole-pairs NS/p
parmot[ite_ratW   ]  = 0.45     # Ratio of tooth width to slot pitch = wT/wST
parmot[ite_ratShft]  = 3/4     # Ratio of inner to outer dia for hollow Shaft
parmot[ite_z      ]  = 3       # Phases

# parmot[ite_kServ  ]  = 1.15    # Rucker2005
# parmot[ite_kServ  ]  = 1.5     # Yoon2016 (Tab.4) Excluding Ground Cylinder and Heat Sink
parmot[ite_kServ  ]  = 1.7     # Yoon2016 (Tab.4) Including Ground Cylinder and Heat Sink

parmot[ite_ratSd ]  = 1/50    # Slot depression to height ratio hsd/hsd
parmot[ite_kpf   ]  = 0.5     # kpf = syrat = Packing factor
parmot[ite_wSd   ]  = 1e-6    # Slot depression width
parmot[ite_Nshrt ]  = 1       # Number of short pitched slots
parmot[ite_hRS   ]  = 0.002   # [m] Thickness of retaining sleeve for Vtip = 200 m/s
parmot[ite_Tarm  ]  = 273.16 + 80 # Armature Temperature in K
parmot[ite_kst   ]  = 0.98    # Stacking factor for silicon electrical steel

parmot[ite_p]     = 8   # pole-pairs
        
parmot[ite_Br  ]  = 1.2 # Magnet remanence (Nd-Fe-B magnets)
parmot[ite_BSat]  = 1.6 # Can go upto 1.9 perhaps
parmot[ite_mu0 ]  = 4π * 1e-7 # Vacuum permeability

parmot[ite_rhoMag  ]     = 7400. # Ne-Fe-B magnet density [kg/m^3]
parmot[ite_rhoCu   ]     = 8900. # Cu density
parmot[ite_rhoFe   ]     = 7700. # Silicon electric Iron density
parmot[ite_rhoSteel]     = 7850. # S275 steel 
parmot[ite_tauMax  ]     = 110.0e6 # Max shaft torsion shear stress S275 Steel
parmot[ite_sigCu   ]     = 1.0e7   # Electrical condutivity of copper at 293K
parmot[ite_thetaCu ]     = 3.81e-3 # Temp coeff of copper (Hanselman)

parmot[ite_psi]     = 0.0 # Power factor angle [rad]

parmot[ite_pb0 ]   = 36.79  # Base factor for power
parmot[ite_Bb0 ]   = 1.0    # Base factor for flux density
parmot[ite_fb0 ]   = 1000.  # Base factor for frequency
parmot[ite_epsb]   = 2.12   # Base factor for flux density exponent
parmot[ite_epsf]   = 1.68   # Base factor for frequency exponent


# ------------- GENERATOR ----------------
# Parameters for the Turbo-electric propulsion system
pargen = zeros(Float64, ite_total)
# PMSM parameters
pargen[ite_ratSM  ]  = 0.86    # Slot to magnet height ratio (hS/hM)
pargen[ite_ratAg  ]  = 1/160   # Air-gap ratio (hAg - hRS)/dRot
pargen[ite_ratM   ]  = 6       # Magnet height ratio hM/hAg
pargen[ite_ratSp  ]  = 36/6    # Number of slots to pole-pairs NS/p
pargen[ite_ratW   ]  = 0.46     # Ratio of tooth width to slot pitch = wT/wST
pargen[ite_ratShft]  = 3/4     # Ratio of inner to outer dia for hollow Shaft
pargen[ite_z      ]  = 3       # Phases

# pargen[ite_kServ  ]  = 1.15    # Rucker2005
# pargen[ite_kServ  ]  = 1.5     # Yoon2016 (Tab.4) Excluding Ground Cylinder and Heat Sink
pargen[ite_kServ  ]  = 1.7     # Yoon2016 (Tab.4) Including Ground Cylinder and Heat Sink

pargen[ite_ratSd ]  = 1/50    # Slot depression to height ratio hsd/hsd
pargen[ite_kpf   ]  = 0.5     # kpf = syrat = Packing factor
pargen[ite_wSd   ]  = 1e-6    # Slot depression width
pargen[ite_Nshrt ]  = 1       # Number of short pitched slots
pargen[ite_hRS   ]  = 0.002   # [m] Thickness of retaining sleeve for Vtip = 200 m/s
pargen[ite_Tarm  ]  = 273.16 + 80 # Armature Temperature in K
pargen[ite_kst   ]  = 0.98    # Stacking factor for silicon electrical steel

pargen[ite_p]     = 8   # pole-pairs
        
pargen[ite_Br  ]  = 1.2 # Magnet remanence (Nd-Fe-B magnets)
pargen[ite_BSat]  = 1.6 # Can go upto 1.9 perhaps
pargen[ite_mu0 ]  = 4π * 1e-7 # Vacuum permeability

pargen[ite_rhoMag  ]     = 7400. # Ne-Fe-B magnet density [kg/m^3]
pargen[ite_rhoCu   ]     = 8900. # Cu density
pargen[ite_rhoFe   ]     = 7700. # Silicon electric Iron density
pargen[ite_rhoSteel]     = 7850. # S275 steel 
pargen[ite_tauMax  ]     = 110.0e6 # Max shaft torsion shear stress S275 Steel
pargen[ite_sigCu   ]     = 1.0e7   # Electrical condutivity of copper at 293K
pargen[ite_thetaCu ]     = 3.81e-3 # Temp coeff of copper (Hanselman)

pargen[ite_psi]     = 0.0 # Power factor angle [rad]

pargen[ite_pb0 ]   = 36.79  # Base factor for power
pargen[ite_Bb0 ]   = 1.0    # Base factor for flux density
pargen[ite_fb0 ]   = 1000.  # Base factor for frequency
pargen[ite_epsb]   = 2.12   # Base factor for flux density exponent
pargen[ite_epsf]   = 1.68   # Base factor for frequency exponent


