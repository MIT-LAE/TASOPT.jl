function powersizing!(ac, engine_type, ip)
    pare = view(ac.pared, :, ip)
    
    Pdes = ac.pared[iePfanmax]
    u_LT = engine.PEMFC_inputs()

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
    α_guess = pare[iealphawat]

    n_cells, A_cell, Q = engine.PEMsize(Pdes, Vdes, u_LT, α_guess)

    pare[ieQheat] = Q
    ac.pared[iencells, :] .= n_cells
    ac.pared[ieAcell, :] .= A_cell
end

function poweroper!(ac, engine_type, ip)
    pare = view(ac.pared, :, ip)
    
    u_LT = engine.PEMFC_inputs()

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
    j_guess = pare[iejdens]
    α_guess = pare[iealphawat]

    P = pare[iePfan]
    mfuel, V_stack, Q, j, α = engine.PEMoper(P, n_cells, A_cell, u_LT, j_guess, α_guess)

    pare[iemfuel] = mfuel * ac.parg[igneng]
    pare[ieVstack] = V_stack
    pare[ieQheat] = Q
    pare[iejdens] = j
    pare[iealphawat] = α
  
end