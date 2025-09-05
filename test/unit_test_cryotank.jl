using StaticArrays

@testset "Cryogenic tank" begin
    p_atm = TASOPT.p_atm
    Δp = 0.001 * p_atm
    @testset "Thermodynamic properties" begin
        #Hydrogen
        species = "H2"

        #Data from NIST for saturated vapor
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [18.232, 20.369, 22.965, 27.314]
        ρs = [0.71653, 1.3322, 2.5133, 6.1715]
        us = [364.55, 372.65, 379.15, 378.1] * 1e3
        hs = [435.25, 448.71, 459.78, 460.2] * 1e3
        
        @testset "Hydrogen gas" begin
            for (i,p) in enumerate(ps)
                Tsat, ρ, ρ_p, h, u, u_p =  TASOPT.CryoTank.gas_properties(species, p)
                @test Tsat ≈ Tsats[i] atol = 0.5
                @test ρ ≈ ρs[i] atol = 1e-1
                @test u ≈ us[i] rtol = 1e-2
                @test h ≈ hs[i] rtol = 1e-2

                #Check derivative by finite differences
                _, ρ2, _, _, u2, _ =  TASOPT.CryoTank.gas_properties(species, p + Δp)
                ρ_pFD = (ρ2 - ρ) / Δp
                u_pFD = (u2 - u) / Δp
                @test ρ_p ≈ ρ_pFD rtol = 1e-2
                @test u_p ≈ u_pFD rtol = 1e-2
            end
        end

        #Data from NIST for saturated liquid
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [18.232, 20.369, 22.965, 27.314]
        ρs = [73.142, 70.848, 67.639, 60.711]
        us = [-20.876, -1.4302, 25.31, 80.495] * 1e3
        hs = [-20.183, 0, 28.306, 88.84] * 1e3

        @testset "Hydrogen liquid" begin
            for (i,p) in enumerate(ps)
                Tsat, ρ, ρ_p, h, u, u_p =  TASOPT.CryoTank.liquid_properties(species, p)
                @test Tsat ≈ Tsats[i] atol = 0.5
                @test ρ ≈ ρs[i] rtol = 1e-2
                @test u ≈ us[i] atol = 2e3
                @test h ≈ hs[i] atol = 2e3

                _, ρ2, _, _, u2, _ =  TASOPT.CryoTank.liquid_properties(species, p + Δp)
                ρ_pFD = (ρ2 - ρ) / Δp
                u_pFD = (u2 - u) / Δp
                @test ρ_p ≈ ρ_pFD rtol = 1e-2
                @test u_p ≈ u_pFD rtol = 1e-2
            end
        end

        #Methane
        species = "CH4"

        #Data from NIST for saturated vapor
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [103.87, 111.67, 120.81, 135.59]
        ρs = [0.96233, 1.8164, 3.4382, 8.1024]
        us = [444.65, 455.05, 466.28, 481.52] * 1e3
        hs = [497.3, 510.83, 525.22, 544.05] * 1e3

        @testset "Methane gas" begin
            for (i,p) in enumerate(ps)
                Tsat, ρ, ρ_p, h, u, u_p =  TASOPT.CryoTank.gas_properties(species, p)
                @test Tsat ≈ Tsats[i] atol = 0.6
                @test ρ ≈ ρs[i] atol = 1e-1
                @test u ≈ us[i] rtol = 1e-2
                @test h ≈ hs[i] rtol = 1e-2

                _, ρ2, _, _, u2, _ =  TASOPT.CryoTank.gas_properties(species, p + Δp)
                ρ_pFD = (ρ2 - ρ) / Δp
                u_pFD = (u2 - u) / Δp
                @test ρ_p ≈ ρ_pFD rtol = 1e-2
                @test u_p ≈ u_pFD rtol = 1e-2
            end
        end

        #Data from NIST for saturated liquid
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [103.87, 111.67, 120.81, 135.59]
        ρs = [433.5, 422.36, 408.66, 384.63]
        us = [-27.112, 0, 31.803, 85.119] * 1e3
        hs = [-26.995, 0, 32.299, 86.436] * 1e3

        @testset "Methane liquid" begin
            for (i,p) in enumerate(ps)
                Tsat, ρ, ρ_p, h, u, u_p =  TASOPT.CryoTank.liquid_properties(species, p)
                @test Tsat ≈ Tsats[i] atol = 0.6
                @test ρ ≈ ρs[i] rtol = 1e-2
                @test u ≈ us[i] atol = 2e3
                @test h ≈ hs[i] atol = 2e3

                _, ρ2, _, _, u2, _ =  TASOPT.CryoTank.liquid_properties(species, p + Δp)
                ρ_pFD = (ρ2 - ρ) / Δp
                u_pFD = (u2 - u) / Δp
                @test ρ_p ≈ ρ_pFD rtol = 1e-2
                @test u_p ≈ u_pFD rtol = 1e-2
            end
        end
    end

    @testset "Saturated mixtures" begin
        species = "H2"
        
        p = p_atm
        β = 0.5

        @testset "Saturated properties" begin
            Tsat, ρ, ρ_p, h, u, u_p = TASOPT.CryoTank.gas_properties(species, p)
            gas = TASOPT.CryoTank.SaturatedGas(species, p)

            @test gas.p ≈ p
            @test gas.T ≈ Tsat
            @test gas.ρ ≈ ρ
            @test gas.u ≈ u
            @test gas.h ≈ h
            @test gas.ρ_p ≈ ρ_p
            @test gas.u_p ≈ u_p

            Tsat, ρ, ρ_p, h, u, u_p = TASOPT.CryoTank.liquid_properties(species, p)
            liq = TASOPT.CryoTank.SaturatedLiquid(species, p)

            @test liq.p ≈ p
            @test liq.T ≈ Tsat
            @test liq.ρ ≈ ρ
            @test liq.u ≈ u
            @test liq.h ≈ h
            @test liq.ρ_p ≈ ρ_p
            @test liq.u_p ≈ u_p

            mix = TASOPT.CryoTank.SaturatedMixture(species, p, β)

            x = 1 / (1 + (liq.ρ/gas.ρ) * (β / (1 - β)))
            ρ = 1 / (x / gas.ρ + (1-x) / liq.ρ)
            u = x * gas.u + (1-x) * liq.u
            h = x * gas.h + (1-x) * liq.h

            x_p = (-x/gas.ρ^2 * gas.ρ_p - (1-x)/liq.ρ^2 * liq.ρ_p) / (1/liq.ρ - 1/gas.ρ)
            u_p = x * gas.u_p + (1 - x) * liq.u_p + x_p * (gas.u - liq.u)
            ϕ = 1 / (u_p * ρ)
            hvap = gas.h - liq.h
            ρ_star = gas.ρ / (liq.ρ - gas.ρ)

            @test mix.p ≈ p
            @test mix.T ≈ Tsat
            @test mix.ρ ≈ ρ
            @test mix.x ≈ x
            @test mix.h ≈ h
            @test mix.u ≈ u
            @test mix.β ≈ β
            @test mix.u_p ≈ u_p
            @test mix.ϕ ≈ ϕ
            @test mix.hvap ≈ hvap
            @test mix.ρ_star ≈ ρ_star
        end
        @testset "Mixture functions" begin     

            mix = TASOPT.CryoTank.SaturatedMixture(species, p, β)

            #Test update to mixture by comparing against new object
            p2 = 2 * p_atm 
            β2 = 0.3
            TASOPT.CryoTank.update_pβ!(mix, p2, β2)

            mix2 = TASOPT.CryoTank.SaturatedMixture(species, p2, β2)

            @test mix2.p ≈ mix.p
            @test mix2.T ≈ mix.T
            @test mix2.ρ ≈ mix.ρ
            @test mix2.x ≈ mix.x
            @test mix2.h ≈ mix.h
            @test mix2.u ≈ mix.u
            @test mix2.β ≈ mix.β
            @test mix2.u_p ≈ mix.u_p
            @test mix2.ϕ ≈ mix.ϕ
            @test mix2.hvap ≈ mix.hvap
            @test mix2.ρ_star ≈ mix.ρ_star

            #Check that density is the same during pressure change with convert_β_same_ρ
            mix = TASOPT.CryoTank.SaturatedMixture(species, p, β)
            ρorig = mix.ρ

            β3 = TASOPT.CryoTank.convert_β_same_ρ(species, p2, p, β)
            TASOPT.CryoTank.update_pβ!(mix, p2, β3)
            @test mix.ρ ≈ ρorig

        end
    end

    species = "h2"
    p = 2.0 * p_atm
    β = 0.5 
    mix = TASOPT.CryoTank.SaturatedMixture("h2", p, β)

    Q = 2e4
    mdot = 1.0
    W = 0.0
    mdot_vent = 0.0

    xout = 0.0 #Taking liquid out
    xvent = 1.0 #Venting gas
    α = 2.0
    V = 100.0
    pmax = 5*p_atm

    @testset "Model derivatives" begin
        #Check function for dp/dt
        dp_dt = TASOPT.CryoTank.dpdt(mix, Q, W, mdot, xout, mdot_vent, xvent, V, α)
        dp_dt_validate = α * mix.ϕ / V * (Q + W - mdot * mix.hvap * (xout + mix.ρ_star) - mdot_vent * mix.hvap * (xvent + mix.ρ_star))
        dp_dt_regression = 6.995703982056127

        @test dp_dt ≈ dp_dt_validate
        @test dp_dt ≈ dp_dt_regression

        #Check function for dβ/dt
        dβ_dt = TASOPT.CryoTank.dβdt(mix, dp_dt_regression, mdot + mdot_vent, V)
        dβ_dt_validate = (-(mdot + mdot_vent) / V - dp_dt_validate * (β * mix.liquid.ρ_p + (1-β)*mix.gas.ρ_p)) / (mix.liquid.ρ - mix.gas.ρ)
        dβ_dt_regression = -0.00015294953663597238

        @test dβ_dt ≈ dβ_dt_validate
        @test dβ_dt ≈ dβ_dt_regression

        #Check function for venting mass flow rate
        mdot_vent_out = TASOPT.CryoTank.venting_mass_flow(mix, Q, W, mdot, xout, xvent) #This venting should give zero dp/dt
        vent_validate = (Q - mdot * mix.hvap * (xout + mix.ρ_star)) / (mix.hvap * (xvent + mix.ρ_star))
        vent_regression = 0.007274844258979718

        @test mdot_vent_out ≈ vent_validate
        @test mdot_vent_out ≈ vent_regression

        dpdt_venting = TASOPT.CryoTank.dpdt(mix, Q, W, mdot, xout, mdot_vent_out, xvent, V, α) #Check pressure derivative with venting; it should be 0
        dpdt_venting_sanity = 0.0

        @test dpdt_venting ≈ dpdt_venting_sanity

        #Check function for boiloff rate
        mdot_liq = (1 - xout) * mdot + (1 - xvent) * mdot_vent
        mdotboil = TASOPT.CryoTank.mdot_boiloff(mix, dβ_dt_regression, dp_dt_regression, mdot_liq, V)
        mdotboil_validate = - (dβ_dt_regression * V * mix.liquid.ρ + β * V * mix.liquid.ρ_p * dp_dt_regression + mdot)
        mdotboil_regression = 0.04266684331352266

        @test mdotboil ≈ mdotboil_validate
        @test mdotboil ≈ mdotboil_regression

        #Check function that computes all derivatives for inputs
        Q_calc(t) = Q
        W_calc(t) = W
        mdot_calc(t) = mdot

        u = TASOPT.CryoTank.tank_inputs(Q_calc, W_calc, mdot_calc)
        params = TASOPT.CryoTank.tank_params(mix, V, pmax, xout, xvent, α)

        y = @SVector [p, β, V * mix.ρ, 0.0, 0.0, 0.0]
        t = 0.0

        dydt = TASOPT.CryoTank.TankDerivatives(t, y, u, params)
        @test dydt[1] ≈ dp_dt_regression
        @test dydt[2] ≈ dβ_dt_regression
        @test dydt[3] ≈ -mdot
        @test dydt[4] ≈ mdot
        @test dydt[5] ≈ 0.0
        @test dydt[6] ≈ mdotboil_regression
    end
end