"""
Thermodynamic properties of one saturated phase (gas or liquid) at a given pressure.

Returned by [`gas_properties`](@ref) and [`liquid_properties`](@ref).

"""
struct SaturatedPhaseProps
    Tsat::Float64   # saturation temperature [K]
    ρ::Float64      # density [kg/m³]
    ρ_p::Float64    # ∂ρ/∂p [kg/(m³·Pa)]
    h::Float64      # specific enthalpy [J/kg]
    u::Float64      # specific internal energy [J/kg]
    u_p::Float64    # ∂u/∂p [J/(kg·Pa)]
end

"""
    gas_properties(species::String, p::Float64)

This function returns the thermodynamic properties of a saturated vapor.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `species::String`: Species name
    - `p::Float64`: pressure (Pa)
    
    **Output:**
    - `props::SaturatedPhaseProps`
        - `props.Tsat::Float64`: temperature (K)
        - `props.ρ::Float64`: density (kg/m^3)
        - `props.ρ_p::Float64`: derivative of density with pressure (kg/m^3/Pa)
        - `props.h::Float64`: specific enthalpy (J/kg)
        - `props.u::Float64`: specific internal energy (J/kg)
        - `props.u_p::Float64`: derivative of internal energy with pressure (J/kg/Pa)
"""
function gas_properties(species::String, p::Float64)
    x = p / p_atm #pressure in atm

    #Fits to NIST data at p increments of 0.1 atm from 0.1 to 10 atm
    if uppercase(species) == "H2" || uppercase(species) == "LH2"
        #Saturation temperature (K)
        Tsat = -3.35772E-04*x^6 + 1.14637E-02*x^5 - 1.54871E-01*x^4 + 1.05803E+00*x^3 - 3.93770E+00*x^2 + 9.09205E+00*x + 1.42913E+01
        
        #Gas density (kg/m^3)
        ρ = 5.61446E-03*x^3 - 4.19141E-02*x^2 + 1.28761E+00*x + 6.91405E-02
        #Density-pressure derivative (kg/m^3/atm)
        ρ_p = 3*5.61446E-03*x^2 - 2*4.19141E-02*x + 1.28761E+00
        #Gas enthalpy (kJ/kg)
        h = -3.31730E-03*x^6 + 1.12617E-01*x^5 - 1.51499E+00*x^4 + 1.02766E+01*x^3 - 3.79360E+01*x^2 + 7.33546E+01*x + 4.04343E+02
        #Gas internal energy (kJ/kg)
        u = -2.01181E-03*x^6 + 6.82454E-02*x^5 - 9.18066E-01*x^4 + 6.22852E+00*x^3 - 2.30505E+01*x^2 + 4.44191E+01*x + 3.45858E+02
        #Internal energy-pressure derivative (kJ/kg/atm)
        u_p = -6*2.01181E-03*x^5 + 5*6.82454E-02*x^4 - 4*9.18066E-01*x^3 + 3*6.22852E+00*x^2 - 2*2.30505E+01*x + 4.44191E+01
    elseif uppercase(species) == "CH4"
        Tsat = -9.38765E-04*x^6 + 3.25637E-02*x^5 - 4.49376E-01*x^4 + 3.16285E+00*x^3 - 1.22945E+01*x^2 + 3.00925E+01*x + 9.09645E+01
        
        ρ = 1.66091E-03*x^3 - 2.82413E-02*x^2 + 1.69557E+00*x + 1.33041E-01
        ρ_p = 3*1.66091E-03*x^2 - 2*2.82413E-02*x + 1.69557E+00

        h = -1.88367E-03*x^6 + 6.52861E-02*x^5 - 8.99795E-01*x^4 + 6.31818E+00*x^3 - 2.44324E+01*x^2 + 5.58852E+01*x + 4.73575E+02
        u = -1.41615E-03*x^6 + 4.90837E-02*x^5 - 6.76548E-01*x^4 + 4.75172E+00*x^3 - 1.83908E+01*x^2 + 4.24831E+01*x + 4.26591E+02
        u_p = -6*1.41615E-03*x^5 + 5*4.90837E-02*x^4 - 4*6.76548E-01*x^3 + 3*4.75172E+00*x^2 - 2*1.83908E+01*x + 4.24831E+01
    elseif uppercase(species) == "C2H6" || uppercase(species) == "ETHANE"
      # webscraper used to quiery 0.1:0.1:10 atm data over a large temperature range to extract saturation conditions for each pressure.
      Tsat = -2.1077168641e-03*x^6 + 7.1830433565e-02*x^5 - 9.6763721978e-01*x^4 + 6.5772804573e+00*x^3 - 2.4201495129e+01*x^2 + 5.3516481777e+01*x + 1.4952928502e+02
      ρ =  1.8394749479e-03*x^3 -3.0510482305e-02*x^2 + 1.9386715191e+00*x + 1.2396908160e-01
      ρ_p = 5.5184248438e-03*x^2 -6.1020964611e-02*x + 1.9386715191e+00
      h =  -2.6664579591e-03*x^6 + 9.0862961217e-02*x^5 -1.2236977065e+00*x^4 + 8.3113102681e+00*x^3 - 3.0514489247e+01*x^2 + 6.5056279260e+01*x + 4.4762676350e+02
      u = -2.1007890231e-03*x^6 + 7.1604017142e-02*x^5 -9.6468295999e-01*x^4 + 6.5564013494e+00*x^3 -2.4110014447e+01*x^2 + 5.2111921806e+01*x + 4.0637144562e+02
      u_p = -1.2604734139e-02*x^5 + 3.5802008571e-01*x^4 -3.8587318399e+00*x^3 + 1.9669204048e+01*x^2 -4.8220028894e+01*x + 5.2111921806e+01
    elseif uppercase(species) == "C2H4" || uppercase(species) == "ETHYLENE"
      Tsat = -1.9243300335e-03*x^6 + 6.5582429082e-02*x^5 - 8.8350874668e-01*x^4 + 6.0059460047e+00*x^3 - 2.2103859000e+01*x^2 + 4.8933219463e+01*x + 1.3732204674e+02
      ρ =  1.8819050342e-03*x^3 - 3.2021905627e-02*x^2 + 1.9701927074e+00*x + 1.2668025371e-01
      ρ_p = 5.6457151026e-03*x^2 -6.4043811253e-02*x + 1.9701927074e+00
      h = -2.2215897848e-03*x^6 + 7.5666036135e-02*x^5 -1.0183183370e+00*x^4 + 6.9088678919e+00*x^3 - 2.5310595095e+01*x^2 + 5.3063591688e+01*x + 4.4864975680e+02
      u = -1.6708380580e-03*x^6 + 5.6912401878e-02*x^5 -7.6603565644e-01*x^4 + 5.1986555487e+00*x^3 - 1.9061036362e+01*x^2 + 4.0365668953e+01*x + 4.0805168296e+02
      u_p = -1.0025028348e-02*x^5 + 2.8456200939e-01*x^4 -3.0641426258e+00*x^3 + 1.5595966646e+01*x^2 -3.8122072723e+01*x + 4.0365668953e+01
    end

    h = h * 1e3 #J/kg
    u = u * 1e3 #J/kg
    ρ_p = ρ_p / p_atm #kg/m^3/Pa
    u_p = u_p * 1e3 / p_atm #J/kg/Pa
    return SaturatedPhaseProps(Tsat, ρ, ρ_p, h, u, u_p)
end

"""
    liquid_properties(species::String, p::Float64)

This function returns the thermodynamic properties of a saturated liquid.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `species::String`: Species name
    - `p::Float64`: pressure (Pa)
    
    **Output:**
    - `props::SaturatedPhaseProps`
        - `props.Tsat::Float64`: temperature (K)
        - `props.ρ::Float64`: density (kg/m^3)
        - `props.ρ_p::Float64`: derivative of density with pressure (kg/m^3/Pa)
        - `props.h::Float64`: specific enthalpy (J/kg)
        - `props.u::Float64`: specific internal energy (J/kg)
        - `props.u_p::Float64`: derivative of internal energy with pressure (J/kg/Pa)
"""
function liquid_properties(species::String, p::Float64)
    x = p / p_atm #pressure in atm

    #Fits to NIST data at p increments of 0.1 atm from 0.1 to 10 atm
    if uppercase(species) == "H2" || uppercase(species) == "LH2"
        #Saturation temperature (K)
        Tsat = -3.35772E-04*x^6 + 1.14637E-02*x^5 - 1.54871E-01*x^4 + 1.05803E+00*x^3 - 3.93770E+00*x^2 + 9.09205E+00*x + 1.42913E+01
        
        #Liquid density (kg/m^3)
        ρ = 2.58354E-04*x^6 - 8.90930E-03*x^5 + 1.21248E-01*x^4 - 8.37637E-01*x^3 + 3.14780E+00*x^2 - 8.42843E+00*x + 7.68629E+01
        #Density-pressure derivative (kg/m^3/atm)
        ρ_p = 6*2.58354E-04*x^5 - 5*8.90930E-03*x^4 + 4*1.21248E-01*x^3 - 3*8.37637E-01*x^2 + 2*3.14780E+00*x - 8.42843E+00
        #Liquid enthalpy (kJ/kg)
        h = -2.258517E-03*x^6 + 7.767483E-02*x^5 - 1.055739E+00*x^4 + 7.280899E+00*x^3 - 2.739412E+01*x^2 + 7.379258E+01*x - 5.277765E+01
        #Liquid internal energy (kJ/kg)
        u = -2.268218E-03*x^6 + 7.793319E-02*x^5 - 1.058739E+00*x^4 + 7.297785E+00*x^3 - 2.749429E+01*x^2 + 7.244553E+01*x - 5.277403E+01
        #Internal energy-pressure derivative (kJ/kg/atm)
        u_p = -6*2.268218E-03*x^5 + 5*7.793319E-02*x^4 - 4*1.058739E+00*x^3 + 3*7.297785E+00*x^2 - 2*2.749429E+01*x + 7.244553E+01
    elseif uppercase(species) == "CH4" || uppercase(species) == "LCH4"
        Tsat = -9.38765E-04*x^6 + 3.25637E-02*x^5 - 4.49376E-01*x^4 + 3.16285E+00*x^3 - 1.22945E+01*x^2 + 3.00925E+01*x + 9.09645E+01
        
        ρ = 1.23474E-03*x^6 - 4.28482E-02*x^5 + 5.91776E-01*x^4 - 4.17205E+00*x^3 + 1.62724E+01*x^2 - 4.14826E+01*x + 4.51398E+02
        ρ_p = 6*1.23474E-03*x^5 - 5*4.28482E-02*x^4 + 4*5.91776E-01*x^3 - 3*4.17205E+00*x^2 + 2*1.62724E+01*x - 4.14826E+01

        h = -3.14958E-03*x^6 + 1.09246E-01*x^5 - 1.50774E+00*x^4 + 1.06164E+01*x^3 - 4.13089E+01*x^2 + 1.02732E+02*x - 7.11737E+01
        u = -3.14796E-03*x^6 + 1.09204E-01*x^5 - 1.50735E+00*x^4 + 1.06151E+01*x^3 - 4.13130E+01*x^2 + 1.02493E+02*x - 7.11704E+01
        u_p = -6*3.14796E-03*x^5 + 5*1.09204E-01*x^4 - 4*1.50735E+00*x^3 + 3*1.06151E+01*x^2 - 2*4.13130E+01*x + 1.02493E+02
    elseif uppercase(species) == "C2H6" || uppercase(species) == "ETHANE"
      # webscraper used to quiery 0.1:0.1:10 atm data over a large temperature range to extract saturation conditions for each pressure.
      Tsat = -2.1077168641e-03*x^6 + 7.1830433565e-02*x^5 - 9.6763721978e-01*x^4 + 6.5772804573e+00*x^3 - 2.4201495129e+01*x^2 + 5.3516481777e+01*x + 1.4952928502e+02
      ρ =  -1.6478244926e-01*x^3 + 3.3782835909e+00*x^2 -2.8594155617e+01*x + 5.7276516047e+02
      ρ_p = -4.9434734779e-01*x^2 + 6.7565671819e+00*x -2.8594155617e+01
      h =  -4.8966018967e-03*x^6 + 1.6691888560e-01*x^5 -2.2495165917e+00*x^4 + 1.5301968437e+01*x^3 -5.6392115502e+01*x^2 + 1.2685386096e+02*x -8.3784297518e+01
      u = -4.8969430560e-03*x^6 + 1.6693104084e-01*x^5 -2.2496930602e+00*x^4 + 1.5303343081e+01*x^3 -5.6400876660e+01*x^2 + 1.2667396286e+02*x -8.3783196657e+01
      u_p = -2.9381658336e-02*x^5 + 8.3465520419e-01*x^4 -8.9987722409e+00*x^3 + 4.5910029243e+01*x^2 -1.1280175332e+02*x + 1.2667396286e+02
    elseif uppercase(species) == "C2H4" || uppercase(species) == "ETHYLENE"
      Tsat = -1.9243300335e-03*x^6 + 6.5582429082e-02*x^5 - 8.8350874668e-01*x^4 + 6.0059460047e+00*x^3 - 2.2103859000e+01*x^2 + 4.8933219463e+01*x + 1.3732204674e+02
      ρ =  -1.7307188938e-01*x^3 + 3.5481537383e+00*x^2 -2.9926640616e+01*x + 5.9794358556e+02
      ρ_p = -5.1921566814e-01*x^2 + 7.0963074766e+00*x -2.9926640616e+01
      h =  -4.6194892733e-03*x^6 + 1.5740220241e-01*x^5 -2.1197821348e+00*x^4 + 1.4401563844e+01*x^3 -5.2925354431e+01*x^2 + 1.1763812134e+02*x -7.7243253370e+01
      u = -4.6198182236e-03*x^6 + 1.5741392408e-01*x^5 -2.1199523083e+00*x^4 + 1.4402890184e+01*x^3 -5.2933768344e+01*x^2 + 1.1746577927e+02*x -7.7242191852e+01
      u_p = -2.7718909341e-02*x^5 + 7.8706962038e-01*x^4 -8.4798092330e+00*x^3 + 4.3208670553e+01*x^2 -1.0586753669e+02*x + 1.1746577927e+02
    end

    h = h * 1e3 #J/kg
    u = u * 1e3 #J/kg
    ρ_p = ρ_p / p_atm #kg/m^3/Pa
    u_p = u_p * 1e3 / p_atm #J/kg/Pa

    return SaturatedPhaseProps(Tsat, ρ, ρ_p, h, u, u_p)
end
