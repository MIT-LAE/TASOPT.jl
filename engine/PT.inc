# Parameters for the PowerTrain
parpt = zeros(Union{Int64, Float64}, ipt_total)

parpt[ipt_nfan   ] = 6
parpt[ipt_ngen   ] = 4
parpt[ipt_nTshaft] = 2
parpt[ipt_Fnsplit] = 1/2

#Turboshaft
parpt[ipt_pifan]   = 1.685
parpt[ipt_piLPC]   = 1.935
parpt[ipt_piHPC]   = 9.369
parpt[ipt_Tt41 ]   = 2750.0

#PCEC
parpt[ipt_cpsi  ]  = 900.0   # Catalyst cells per square inch (CPSI)
parpt[ipt_wcat  ]  = 2.0     # Catalyst wall thickness in [mil]
parpt[ipt_lcat  ]  = 0.02    # Catalyst length
parpt[ipt_deNOx ]  = 0.98    # DeNOx target at cruise

#Generator
parpt[ipt_ARgen]       = 0.6
parpt[ipt_sigAgGen]    = 45e3
parpt[ipt_ratSplitGen] = 0.82

#Motors
parpt[ipt_ARmot]       = 0.6
parpt[ipt_sigAgMot]    = 30e3
parpt[ipt_ratSplitMot] = 0.7

#Cables
parpt[ipt_Vcable] = 540 # Volts
parpt[ipt_sigcon  ] = 5.96e7 # Conductivity of Copper @ 293K [S/m]
parpt[ipt_alphacon] = 0.00429
parpt[ipt_rhocon]   = 8960.0 # Density of copper
parpt[ipt_Jmax]     = 5.0e6  #A/m²
parpt[ipt_kpf]      = 0.91 # Packing factor of Litz cables from Dowdle 2018
parpt[ipt_rhoins]   = 1700 #[kg/m³] polyamide insulation #Dowdle 2018
parpt[ipt_Emax]     = 1.0e7 #V/m Polyamide dielectric strength Dowdle 2018

# Weights are also stored here
