include("wingpo.jl")
include("atmos.jl")
include("tankWmech.jl")
include("tankWthermal.jl")

print(varinfo())

(gee, Nland,
      Wfix, Wpay, Wpadd, Wseat, Wapu, Weng,
      fstring, fframe, ffadd,
      deltap,
      Wpwindow, Wppinsul, Wppfloor,
      Whtail, Wvtail,
      rMh, rMv,
      Lhmax, Lvmax,
      bv, lambdav, nvtail,
      Rfuse, dRfuse, wfb, nfweb,
      lambdac, xnose, xshell1, xshell2, xconend,
      xhtail, xvtail,
      xwing, xwbox,
      cbox, xfix, xapu, xeng,
      hfloor, sigskin, sigbend,
      rhoskin, rhobend,
      Eskin, Ebend, Gskin, rhoFuel, m_airplane, R, lcv, eta, LD, v_cruise, hconvgas,
      h_LH2, Tfuel, Tair, r_tank, h_v, t_cond, r_gas, k, hconvair, time_flight, Wfuel, lshell, threshold_percent, mode) = (9.81, 6.0,
      0.13E+05,  0.46E+06,  0.16E+06,  0.46E+05,  0.16E+05,   0.0,
      0.34,      0.24,      0.20,
      0.56E+05,
      0.14E+03,   40,       60,      0.26E+05,  0.22E+05,
      0.40,      0.70,      0.17E+07,  0.15E+07,
      12,      0.25,       1.0,       3.1,
      0.0,       0.0,       1.0,      0.30,
      2.4,       12.,       62.,       72,
      70.,       67.,       40.,       34,       5.7,       3.0,       71.,       31.,
      0.20,      0.10E+09,  0.21E+09,  0.27E+04,
      0.27E+04,  0.69E+11,  0.69E+11,  0.27E+11, 67.3, 33000, 5000000, 80000000, 0.5, 18, 280, 25,
      50, 21, 173, 3, 80000000, [0.1, 0.1, 0.1], 0.1, [205.0, 205.0, 205.0], 100, 3600, 5000, 2, 5, 2)

rho_insul = [100, 100, 100]
t=[.1 .1 .1]
#rhoFuel from hydrogen tank designv3 article
thickness_insul = sum(t_cond)
m_boiloff = tankWthermal(lshell, hconvgas, h_LH2, Tfuel, Tair, r_tank,
                      h_e, t_cond, k, hconvair, time_flight)
result = tankWmech(gee, rhoFuel,
                      fstring, ffadd, deltap,
                      Rfuse, dRfuse, wfb, nfweb,
                      sigskin, rho_insul, rhoskin,
                      Wfuel, m_boiloff, thickness_insul, t_cond)

Wtank = result[1]
#Wfuel = result[2]

#Wtank = Wtank + m_boiloff*gee + Wfuel #weight of tank including fuel

#print(result)
#print("atmos(0)")
#print("atmos vars = ", atmos(0))
#a=atmos(0)
#print(a)
