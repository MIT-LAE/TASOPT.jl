include("wingpo.jl")
include("fusew.jl")

#print(varinfo())

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

include("blsys.jl")

(simi, lami, wake, direct) = (false, false, true, true)
Mach = 0.83999999999999997 
uinv = 0.99923266171096192
hksep = 2.8999999999999999

x = 113.19176121887341
b = 4.2512954905511248E-002 
rn = 1.0000000000000000   
th = 0.36151524580335676 
ds = 0.47411597247856280    
ue = 0.99923266171096192

h = 1.3114688190396659
h_th = -3.6276998944410455
h_ds = 2.7661350706740091 

hk = 1.0255741688623847
hk_th = -3.3602469158702317  
hk_ds = 2.5622011496474624
hk_ue = -0.60482391056544504 

hc = 0.18475056251156852
hc_th = -0.51104501029000049  
hc_ds = 0.38967377864478508 
hc_ue = 0.42195763911809053  

hs = 1.9793306587190489  
hs_th = 2.6598043330637906
hs_ds = -2.0281130759023553 
hs_ue = 0.47920941853514037 

cf = 0.0000000000000000  
cf_th = 0.0000000000000000  
cf_ds = 0.0000000000000000   
cf_ue = 9.7873179306789851E-005

di = 1.1201590943705716E-006
di_th = -4.3053013257921013E-004
di_ds = 3.2828087586136393E-004 
di_ue = 1.6832700785353640E-005

xm = 109.41631855980170 
bm = 4.2512954905511248E-002
rnm = 1.0000000000000000   
thm = 0.36159007642163560 
dsm = 0.47420811585137385
uem = 0.99907993608762147 

(hm, hm_thm, hm_dsm) = (1.3114522484669959, -3.6269034544711247, 2.7655627253761956)

(hkm, hkm_thm, hkm_dsm, hkm_uem) = (1.0256511654878535, -3.3595955637851405, 2.5617368590524694, -0.60473072971602160)

(hcm, hcm_thm, hcm_dsm, hcm_uem) = (0.18468380516008143,-0.51075464753136102,0.38945729677035568, 0.42186741128515043)

(hsm, hsm_thm, hsm_dsm, hsm_uem) = (1.9792696450491154,2.6589664547453919, -2.0274998135942890,0.47907891555231130)

(cfm, cfm_thm, cfm_dsm, cfm_uem) = (0.00, 0.00, 0.00, 0.00)

(dim, dim_thm, dim_dsm, dim_uem) = (1.1300768260643012E-006,  -4.3292178714927860E-004, 
                                    3.3010869260038754E-004,   1.6387793785720393E-005)


## aa, bb, rr {FORTRAN OUTPUT SO COLOUMN MAJOR}
#2.7658578153890452        1.3463113092256376        0.0000000000000000 
#2.0797231215108840       -1.0265560356304859        0.0000000000000000
#2.6088376442840238       0.11739903875392768        1.0000000000000000    
#
#-2.7658398168637897       -1.3409512628422848        0.0000000000000000  
#-2.0789019632134287        1.0225066805263412        0.0000000000000000     
#-2.6097422210864201      -0.11717682903140124        0.0000000000000000     
#
#-3.9769711976953351E-008  -1.2241218111693874E-009   0.0000000000000000 

(aa, bb, rr) = blsys(simi,lami,wake,direct, Mach, uinv,hksep, x,b,rn,th,ds,ue,
                      h , h_th, h_ds,
                      hk, hk_th, hk_ds, hk_ue,
                      hc, hc_th, hc_ds, hc_ue,
                      hs, hs_th, hs_ds, hs_ue,
                      cf, cf_th, cf_ds, cf_ue,
                      di, di_th, di_ds, di_ue,
                      xm,bm,rnm,thm,dsm,uem, 
                      hm , hm_thm, hm_dsm,
                      hkm, hkm_thm, hkm_dsm, hkm_uem,
                      hcm, hcm_thm, hcm_dsm, hcm_uem,
                      hsm, hsm_thm, hsm_dsm, hsm_uem,
                      cfm, cfm_thm, cfm_dsm, cfm_uem,
                      dim, dim_thm, dim_dsm, dim_uem)

#2.7658578153890452        1.3463113092256376        0.0000000000000000      
#0.75192698263064928       -2.0388838360580386        0.0000000000000000     
#0.94322912398772940       0.56525093665460657        1.0000000000000000
#
# rr Answer in blsys
#-7.6910181160795242E-009  -8.8941884720010155E-009   0.0000000000000000 

solved_bl_sys = aa\rr

print("\nSolved BL system = ", solved_bl_sys)
print("\nTASOPT FORTRAN output = ", "-7.6910181160795242E-009  -8.8941884720010155E-009   0.0000000000000000" )
