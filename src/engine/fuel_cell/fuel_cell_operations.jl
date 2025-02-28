function size_fuel_cell!(ac, ip, imission)
    pare = view(ac.pare, :, ip, imission)
    fcdata = ac.engine.data
    
    Pdes = fcdata.design_power
    u_LT = engine.PEMFC_inputs()

    u_LT.j = fcdata.current_density[ip, imission]
    u_LT.T = fcdata.FC_temperature[ip, imission]	
    u_LT.p_A = fcdata.FC_pressure[ip, imission]
    u_LT.p_C = fcdata.FC_pressure[ip, imission]
    u_LT.x_H2O_A = pare[iexwanode]
    u_LT.x_H2O_C = pare[iexwcathode]
    u_LT.λ_H2 = pare[ielambdaw]
    u_LT.λ_O2 = pare[ielambdaox]
    u_LT.t_M = pare[ietmembrane]
    u_LT.t_A = pare[ietanode]
    u_LT.t_C = pare[ietcathode]
    u_LT.type = fcdata.type

    Vdes = pare[ieVstack]
    α_guess = pare[iealphawat]

    n_cells, A_cell, Q = engine.PEMsize(Pdes, Vdes, u_LT, α_guess)

    pare[ieQheat] = Q
    ac.pared[iencells, :] .= n_cells
    ac.pared[ieAcell, :] .= A_cell
end

function operate_fuel_cell!(ac, ip)
    pare = view(ac.pared, :, ip)
    
    u_LT = engine.PEMFC_inputs()

    u_LT.T = pare[ieTfc]
    u_LT.p_A = pare[iepfc]
    u_LT.p_C = pare[iepfc]
    u_LT.x_H2O_A = pare[iexwanode]
    u_LT.x_H2O_C = pare[iexwcathode]
    u_LT.λ_H2 = pare[ielambdaw]
    u_LT.λ_O2 = pare[ielambdaox]
    u_LT.t_M = pare[ietmembrane]
    u_LT.t_A = pare[ietanode]
    u_LT.t_C = pare[ietcathode]
    u_LT.type = "HT-PEMFC"

    n_cells = pare[iencells]
    A_cell = pare[ieAcell]
    j_guess = pare[iejdens]
    α_guess = pare[iealphawat]

    P = pare[iePfc]
    mfuel, V_stack, Q, j, α = engine.PEMoper(P, n_cells, A_cell, u_LT, j_guess, α_guess)

    pare[iemfuel] = mfuel * ac.parg[igneng]
    pare[ieVstack] = V_stack
    pare[ieQheat] = Q
    pare[iejdens] = j
    pare[iealphawat] = α
  
end