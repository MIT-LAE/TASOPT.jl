
# Parameters for the Turbo-electric propulsion system
  
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

parte[ite_p]    = 8 # pole-pairs
