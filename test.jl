include("wingpo.jl")

include("fusew.jl")

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
      Eskin, Ebend, Gskin) = (9.81, 6.0,  
      0.13E+05,  0.46E+06,  0.16E+06,  0.46E+05,  0.16E+05,   0.0,
      0.34,      0.24,      0.20,    
      0.56E+05,
      0.14E+03,   40.,       60.,      0.26E+05,  0.22E+05,
      0.40,      0.70,      0.17E+07,  0.15E+07,
      12.,      0.25,       1.0,       3.1,
      0.0,       0.0,       1.0,      0.30,
      2.4,       12.,       62.,       72., 
      70.,       67.,       40.,       34.,       5.7,       3.0,       71.,       31.,
      0.20,      0.10E+09,  0.21E+09,  0.27E+04,
      0.27E+04,  0.69E+11,  0.69E+11,  0.27E+11,)


result = fusew(gee, Nland, 
		    Wfix, Wpay, Wpadd, Wseat, Wapu, Weng, 
		    fstring, fframe, ffadd, 
		    deltap, 
		    Wpwindow, Wppinsul, Wppfloor, Whtail, Wvtail, 
		    rMh, rMv, 
		    Lhmax, Lvmax, bv, lambdav, nvtail, 
		    Rfuse, dRfuse, wfb, nfweb, 
		    lambdac, xnose, xshell1, xshell2, xconend, 
		    xhtail, xvtail, xwing, xwbox, cbox, xfix, xapu, xeng, 
		    hfloor, sigskin, sigbend, 
		    rhoskin, rhobend, 
		    Eskin, Ebend, Gskin) 

print(result)
