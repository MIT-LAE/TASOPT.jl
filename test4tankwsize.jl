include("wingpo.jl")
include("atmos.jl")
include("tankWmech.jl")
include("tankWthermal.jl")
include("tankwsize.jl")

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
      h_LH2, Tfuel, Tair, r_tank, h_e, t, r_gas, k, hconvair) = (9.81, 6.0,
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
      50, 21, 173, 3, 80000000, [0.1, 0.1, 0.1], 0.1, [205.0, 205.0, 205.0], 100)

#rhoFuel from hydrogen tank designv3 article
result = tankWmech(gee, rhoFuel,
                      fstring,fframe,ffadd,deltap,
                      Rfuse,dRfuse,wfb,nfweb,
                      xshell1,xshell2,
                      sigskin,Wppinsul, rhoskin,
                      m_airplane, R, lcv, eta, LD)

Wtank = result[1]

m_boiloff = tankWthermal(gee, rhoFuel, deltap,
                        Rfuse, dRfuse,
                        xshell1, xshell2, v_cruise, hconvgas, h_LH2, Tfuel, Tair, r_tank,
                        h_e, t, r_gas, k, hconvair)

Wtank = Wtank + m_boiloff

#print(result)
#print("atmos(0)")
#print("atmos vars = ", atmos(0))
#a=atmos(0)
#print(a)
