using ForwardDiff
"""
    mcool(ncrowx, Tmrow, Tt3, Tt4, dTstreak, Trrat, efilm, tfilm, StA)

Calculates cooling mass flow requirement.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ncrowx`:    dimension of `Tmrow(.)` and `epsrow(.)` arrays (max number of blade rows)
    - `Tmrow(.)`:  design metal temperature for each blade row
    - `Tt3`:       cooling flow temperature
    - `Tt4`:       hot gas temperature from burner
    - `dTstreak`:  hot-streak temperature increase over Tt4, for first blade row 
    - `Trrat`:     static temperature ratio across each blade row, T4.1 / T4
    - `efilm`:     cooling efficiency = `(Tco-Tci)/(Tmetal-Tci)`
    - `tfilm`:     film effectiveness = `(Tgas-Tfaw)/(Tgas-Tco)` where
      - `Tco` = temperature of cooling air exiting  blade
      - `Tci` = temperature of cooling air entering blade
      - `Tfaw` = film adiabatic wall temperature (for insulated-wall case)
    - StA`:        area-weighted external Stanton number = `St (Asurf/Aflow) cpgas/cpcool`

    **Output:**
    - `ncrow`:      number of blade rows which need cooling
    - `epsrow(.)`:  cooling mass flow ratio for each blade row, `m_c_row/m_air`

"""
function mcool(ncrowx,
      Tmrow, Tt3, Tt4, dTstreak, Trrat,
      efilm, tfilm, StA)

      #---- assume zero cooling mass flow for all blade rows
      epsrow = zeros(ncrowx)
      epsrow_Tt3 = zeros(ncrowx)
      epsrow_Tt4 = zeros(ncrowx)
      epsrow_Trr = zeros(ncrowx)
      for icrow = 1:ncrowx
            epsrow[icrow] = 0.0
            epsrow_Tt3[icrow] = 0.0
            epsrow_Tt4[icrow] = 0.0
            epsrow_Trr[icrow] = 0.0
      end

      #---- step downstream, until cooling is no longer needed
      ncrow = 0
      for icrow = 1:ncrowx
            if (icrow == 1)
                  Tg = Tt4 + dTstreak
                  Tg_Tt4 = 1.0
                  Tg_Trr = 0.0
            else
                  Tg = Tt4 * Trrat^(icrow - 1)
                  Tg_Tt4 = Trrat^(icrow - 1)
                  Tg_Trr = float(icrow - 1) * Tg / Trrat
            end
            #
            #------ cooling effectiveness for this blade row
            theta = (Tg - Tmrow[icrow]) / (Tg - Tt3)
            theta_Tt3 = theta / (Tg - Tt3)
            theta_Tt4 = (1.0 - theta) / (Tg - Tt3) * Tg_Tt4
            theta_Trr = (1.0 - theta) / (Tg - Tt3) * Tg_Trr

            eps0 = StA * (theta * (1.0 - efilm * tfilm) - tfilm * (1.0 - efilm)) /
                   (efilm * (1.0 - theta))
            eps0_theta = StA * (1.0 - efilm * tfilm) /
                         (efilm * (1.0 - theta)) +
                         eps0 / (1.0 - theta)
            eps0_Tt3 = eps0_theta * theta_Tt3
            eps0_Tt4 = eps0_theta * theta_Tt4
            eps0_Trr = eps0_theta * theta_Trr

            if (eps0 < 0.0)
                  return ncrow, epsrow, epsrow_Tt3, epsrow_Tt4, epsrow_Trr
            end
            ncrow = icrow
            epsrow[icrow] = eps0 / (1.0 + eps0)
            epsrow_eps0 = (1.0 - epsrow[icrow]) / (1.0 + eps0)

            epsrow_Tt3[icrow] = epsrow_eps0 * eps0_Tt3
            epsrow_Tt4[icrow] = epsrow_eps0 * eps0_Tt4
            epsrow_Trr[icrow] = epsrow_eps0 * eps0_Trr

      end

      return ncrow, epsrow, epsrow_Tt3, epsrow_Tt4, epsrow_Trr
end # mcool

"""
    Tmcalc(ncrowx, ncrow, Tt3, Tt4, dTstreak, Trrat, efilm, tfilm, StA, epsrow)

Calculates metal temperature for blade row

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
    - `ncrowx`:    dimension of Tmrow(.),epsrow(.) arrays (max number of blade rows)
    - `ncrow`:     number of blade rows which are cooled
    - `epsrow(.)`:  cooling mass flow ratio for each blade row, m_c_row/m_air
    - `Tt3`:       cooling flow temperature
    - `Tt4`:       hot gas temperature from burner
    - `dTstreak`:  hot-streak temperature increase over Tt4, for first blade row 
    - `Trrat`:     static temperature ratio across each blade row, T4.1 / T4
    - `efilm`:     cooling efficiency = (Tco-Tci)/(Tmetal-Tci)
    - `tfilm`:     film effectiveness = (Tgas-Tfaw)/(Tgas-Tco)
                   Tco = temperature of cooling air exiting  blade
                   Tci = temperature of cooling air entering blade
                   Tfaw = film adiabatic wall temperature (for insulated-wall case)
    - `StA`:       area-weighted external Stanton number = St (Asurf/Aflow) cpgas/cpcool
       
      **Output:**
    - `Tmrow(.)`:  design metal temperature for each blade row
"""
function Tmcalc(ncrowx, ncrow,
      Tt3, Tt4, dTstreak, Trrat,
      efilm, tfilm, StA, epsrow)

      # Zygote
      # Tmrow = zeros(ncrowx)
      # buf = Zygote.Buffer(Tmrow, ncrowx)
      # for icrow = 1:ncrow
      #       if (icrow == 1)
      #             Tg = Tt4 + dTstreak
      #       else
      #             Tg = Tt4 * Trrat^(icrow - 1)
      #       end

      #       eps = epsrow[icrow]
      #       theta = (eps * efilm + StA * tfilm * (1.0 - efilm) * (1.0 - eps)) /
      #               (eps * efilm + StA * (1.0 - efilm * tfilm) * (1.0 - eps))

      #       # To make Zygote.jl happy
      #       buf[icrow] = Tg - (Tg - Tt3) * theta
      #       theta = (Tg - buf[icrow]) / (Tg - Tt3)

      # end

      # Tmrow = copy(buf)

      # General
      # Tmrow = zeros(ncrowx)
      # for icrow = 1:ncrow
      #       if (icrow == 1)
      #             Tg = Tt4 + dTstreak
      #       else
      #             Tg = Tt4 * Trrat^(icrow - 1)
      #       end

      #       eps = epsrow[icrow]
      #       theta = (eps * efilm + StA * tfilm * (1.0 - efilm) * (1.0 - eps)) /
      #               (eps * efilm + StA * (1.0 - efilm * tfilm) * (1.0 - eps))

      #       # Original code, Zygote.jl does not like it
      #       Tmrow[icrow] = Tg - (Tg - Tt3) * theta
      #       theta = (Tg - Tmrow[icrow]) / (Tg - Tt3)

      # end

      # ForwardDiff
      if (typeof(Tt3) <: ForwardDiff.Dual || typeof(Tt4) <: ForwardDiff.Dual || typeof(dTstreak) <: ForwardDiff.Dual || typeof(Trrat) <: ForwardDiff.Dual || typeof(efilm) <: ForwardDiff.Dual || typeof(tfilm) <: ForwardDiff.Dual || typeof(StA) <: ForwardDiff.Dual || typeof(epsrow[1]) <: ForwardDiff.Dual)
            T = ForwardDiff.Dual
      else
            T = typeof(Tt3)
      end

      Tmrow = zeros(T, ncrowx)
      for icrow = 1:ncrow
            if (icrow == 1)
                  Tg = Tt4 + dTstreak
            else
                  Tg = Tt4 * Trrat^(icrow - 1)
            end

            eps = epsrow[icrow]
            theta = (eps * efilm + StA * tfilm * (1.0 - efilm) * (1.0 - eps)) /
                    (eps * efilm + StA * (1.0 - efilm * tfilm) * (1.0 - eps))

            # Original code, Zygote.jl does not like it
            Tmrow[icrow] = Tg - (Tg - Tt3) * theta
            theta = (Tg - Tmrow[icrow]) / (Tg - Tt3)

      end


      return Tmrow
end # Tmcalc