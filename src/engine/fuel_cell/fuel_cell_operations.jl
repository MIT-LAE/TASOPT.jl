function size_fuel_cell!(ac, ip, imission)
    fcdata = ac.engine.data
    
    Pdes = fcdata.design_power
    u_LT = engine.PEMFC_inputs()

    u_LT.j = fcdata.current_density[ip, imission]
    u_LT.T = fcdata.FC_temperature[ip, imission]	
    u_LT.p_A = fcdata.FC_pressure[ip, imission]
    u_LT.p_C = fcdata.FC_pressure[ip, imission]
    u_LT.x_H2O_A = fcdata.water_concentration_anode[ip, imission]
    u_LT.x_H2O_C = fcdata.water_concentration_cathode[ip, imission]
    u_LT.λ_H2 = fcdata.λ_H2[ip, imission]
    u_LT.λ_O2 = fcdata.λ_O2[ip, imission]
    u_LT.t_M = fcdata.thickness_membrane
    u_LT.t_A = fcdata.thickness_anode
    u_LT.t_C = fcdata.thickness_cathode
    u_LT.type = fcdata.type

    Vdes = fcdata.design_voltage
    #α_guess = fcdata.α_water[ip, imission]

    #n_cells, A_cell, Q = engine.PEMsize(Pdes, Vdes, u_LT, α_guess)
    n_cells, A_cell, Q = engine.PEMsize(Pdes, Vdes, u_LT)

    fcdata.FC_heat[ip, imission] = Q
    fcdata.number_cells = n_cells
    fcdata.area_cell = A_cell
end

function operate_fuel_cell!(ac, ip, imission)
    pare = view(ac.pare, :, ip, imission)
    fcdata = ac.engine.data
    
    u_LT = engine.PEMFC_inputs()

    u_LT.T = fcdata.FC_temperature[ip, imission]	
    u_LT.p_A = fcdata.FC_pressure[ip, imission]
    u_LT.p_C = fcdata.FC_pressure[ip, imission]
    u_LT.x_H2O_A = fcdata.water_concentration_anode[ip, imission]
    u_LT.x_H2O_C = fcdata.water_concentration_cathode[ip, imission]
    u_LT.λ_H2 = fcdata.λ_H2[ip, imission]
    u_LT.λ_O2 = fcdata.λ_O2[ip, imission]
    u_LT.t_M = fcdata.thickness_membrane
    u_LT.t_A = fcdata.thickness_anode
    u_LT.t_C = fcdata.thickness_cathode
    u_LT.type = fcdata.type

    n_cells = fcdata.number_cells
    A_cell = fcdata.area_cell
    #j_guess = fcdata.current_density[ip, imission]
    #α_guess = fcdata.α_water[ip, imission]

    P = fcdata.FC_power[ip, imission]
    #mfuel, V_stack, Q, j, α = engine.PEMoper(P, n_cells, A_cell, u_LT)
    V_stack, Q = engine.PEMoper(P, n_cells, A_cell, u_LT)

    M_h2 = 2.016e-3 #kg/mol
    F = 96485.3329 # C/mol, Faraday contant
    I = P / V_stack #Current
    mfuel = I * M_h2 / (2*F) #Molar flux of hydrogen gas

    pare[iemfuel] = mfuel * ac.parg[igneng]
    fcdata.stack_voltage[ip, imission] = V_stack
    fcdata.FC_heat[ip, imission] = Q
    fcdata.current_density[ip, imission] = I/A_cell/n_cells
    #fcdata.α_water[ip, imission] = α
  
end