# function powercalc!(ac, engine_type, ip)
# end
engine_type = ""
pare = ac.pared
pare[iejdens,:] .= 1e4
pare[ieTfc,:] .= 353.15
pare[iepanode,:] .= 3e5
pare[iepcathode,:] .= 3e5
pare[iexwanode,:] .= 0.1
pare[iexwcathode,:] .= 0.1
pare[ielambdaw,:] .= 3.0
pare[ielambdaox,:] .= 3.0
pare[ietmembrane,:] .= 100e-6
pare[ietanode,:] .= 250e-6
pare[ietcathode,:] .= 250e-6
pare[ieVstack,:] .= 200

function powersizing!(ac, engine_type, ip)
    pare = view(ac.pared, :, ip)
    
    Pdes = ac.pared[iePfanmax]
    u_LT = TASOPT.engine.PEMFC_inputs()

    u_LT.j = pare[iejdens]
    u_LT.T = pare[ieTfc]
    u_LT.p_A = pare[iepanode]
    u_LT.p_C = pare[iepcathode]
    u_LT.x_H2O_A = pare[iexwanode]
    u_LT.x_H2O_C = pare[iexwcathode]
    u_LT.λ_H2 = pare[ielambdaw]
    u_LT.λ_O2 = pare[ielambdaox]
    u_LT.t_M = pare[ietmembrane]
    u_LT.t_A = pare[ietanode]
    u_LT.t_C = pare[ietcathode]
    u_LT.type = "LT-PEMFC"

    Vdes = pare[ieVstack]

    mfuel, n_cells, A_cell, Q = TASOPT.engine.PEMsize(Pdes, Vdes, u_LT)

    pare[iemfuel] = mfuel
    pare[ieQheat] = Q
    ac.pared[iencells, :] .= n_cells
    ac.pared[ieAcell, :] .= A_cell
end

function poweroper!(ac, engine_type, ip)
    pare = view(ac.pared, :, ip)
    
    u_LT = TASOPT.engine.PEMFC_inputs()

    u_LT.T = pare[ieTfc]
    u_LT.p_A = pare[iepanode]
    u_LT.p_C = pare[iepcathode]
    u_LT.x_H2O_A = pare[iexwanode]
    u_LT.x_H2O_C = pare[iexwcathode]
    u_LT.λ_H2 = pare[ielambdaw]
    u_LT.λ_O2 = pare[ielambdaox]
    u_LT.t_M = pare[ietmembrane]
    u_LT.t_A = pare[ietanode]
    u_LT.t_C = pare[ietcathode]
    u_LT.type = "LT-PEMFC"

    n_cells = pare[iencells]
    A_cell = pare[ieAcell]

    ac.pared[ietcathode, :] .= n_cells
    ac.pared[ieAcell, :] .= A_cell

    P = pare[iePfan]
    mfuel, V_stack, Q = TASOPT.engine.PEMoper(P, n_cells, A_cell, u_LT)

    pare[iemfuel] = mfuel
    pare[ieVstack] = V_stack
    pare[ieQheat] = Q
  
end