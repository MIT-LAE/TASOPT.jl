# Parameters for the PowerTrain
parpt = zeros(Union{Int64, Float64}, ipt_total)

parpt[ipt_nfan   ] = 14
parpt[ipt_ngen   ] = 4
parpt[ipt_nTshaft] = 2

#Turboshaft
parpt[ipt_pifan]   = 1.3
parpt[ipt_piLPC]   = 3.0
parpt[ipt_piHPC]   = 7.0
parpt[ipt_Tt41 ]   = 3200.0

#PCEC
parpt[ipt_cpsi  ]  = 900.0   # Catalyst cells per square inch (CPSI)
parpt[ipt_wcat  ]  = 2.0     # Catalyst wall thickness in [mil]
parpt[ipt_lcat  ]  = 0.02    # Catalyst length
parpt[ipt_deNOx ]  = 0.98    # DeNOx target at cruise

#Generator
parpt[ipt_ARgen]       = 0.8
parpt[ipt_sigAgGen]    = 55e3
parpt[ipt_ratSplitGen] = 0.82

#Motors
parpt[ipt_ARmot]       = 0.7
parpt[ipt_sigAgMot]    = 30e3
parpt[ipt_ratSplitMot] = 0.7


# Weights are also stored here
