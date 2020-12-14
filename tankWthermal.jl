"""
tankWthermal calculates boil-off rate of LH2

Inputs:
-Thermal conductivity array k comprising of k value for each MLI layer
-hconvgas is convective coefficient of insulating purged gas (e.g. N2)
-hconvair is convective coefficient of ambient air
-Thickness array t corresponds to thickness value of each layer in MLI
-h_LH2 corresponds to LH2 convective coefficient
-Tfuel is fuel temperature
-Tair is ambient temperature
-r_tank is tank outer radius
-pi is value of pi
-h_e is heat of combustion of liquid hydrogen
-r_gas is inner radius of gas-purged chamber


Outputs:
- m_boiloff is the boiloff LH2 mass for given mission
"""
function tankWthermal(gee, rhoFuel, deltap,
                      Rfuse, dRfuse,
                      xshell1, xshell2, hconvair, hconvgas, h_LH2, Tfuel, Tair, r_tank, pi,
                      h_e, t, r_gas, k)

#--- effective pressure-vessel length
      lshell = xshell2 - xshell1  #pressure vessel length, can be based on tank volume required to store LH2
      N = size(t) #Number of layers in MLI
      thickness = sum(t) #total thickness of MLI

#--- Heat flux and resistances
      deltaT = Tair - Tfuel
      Rair = 1 / (hconvair * 2 * pi * r_tank * lshell)
      r_inner = r_tank - thickness
      Rgas = 1 / (hconvgas * 2 * pi * (r_inner + r_gas) * lshell)
      R_LH2 = 1 / (h_LH2 * 2 * pi * r_inner * lshell)

      rang = [1, 2, 3]
      R_mli = zeros(N)
      for n in rang
            R_mli[n] = log((r_inner  + t[n])/ (r_inner)) / (2 * pi * lshell * k[n])
            r_inner = r_inner + t[n]
      end

      R_mli = sum(R_mli)
      Req = R_mli + Rair + Rgas + R_LH2

      q = deltaT / Req
      m_boiloff = q / h_e

return  m_boiloff
end
