using Pkg, Dates

println(today())
println("Current location $(pwd())")
using TASOPT
const aerodynamics = TASOPT.aerodynamics
include(TASOPT.__TASOPTindices__)
nmisx = 1

using Profile
using BenchmarkTools

println("Loading aircraft model...")

include(joinpath(TASOPT.__TASOPTroot__, "../test/default_sized.jl"))
ac = load_default_model()
size_aircraft!(ac; printiter = false )

println("\nNotes (from BenchmarkTools Manual):
- The minimum is a robust estimator for the location parameter of the
  time distribution, and should not be considered an outlier
- The median, as a robust measure of central tendency,
  should be relatively unaffected by outliers
- The mean, as a non-robust measure of central tendency,
  will usually be positively skewed by outliers
- The maximum should be considered a primarily noise-driven outlier,
  and can change drastically between benchmark trials.\n\n")

println("Start Benchmarking...")

function benchmark_fuselage_drag()
    println("---------------------------------------")
    println("Fuselage boundary layer calculations")
    println("---------------------------------------")

    #Inputs from Fortran
    xnose   = 0.000   
    xend    = 37.79    
    xblend1 = 6.096    
    xblend2 = 29.56    
    Sfuse = 13.50     
    anose = 1.649     
    btail = 2.000     
    ifclose = 0 
    Mach = 0.80
    nc = 30 
    nbldim = 60
    xbl = zeros(nbldim)
    zbl = zeros(nbldim)
    sbl = zeros(nbldim)
    dybl = zeros(nbldim)
    uinv = zeros(nbldim)

    println("Benchmarking... _axisymm_flow")

    bench = @benchmarkable aerodynamics._axisymm_flow($xnose,$xend,$xblend1,$xblend2,
    $Sfuse, $anose, $btail, $ifclose,
    $Mach, $nc, $nbldim,  $xbl, $zbl, $sbl, $dybl, $uinv) seconds=30 evals=50
    bench_axisymm_flow = run(bench)

    #results
    nbl =           47 ;iblte =          31 ;
    # _axisymm_BL benchmarks
    # Load required inputs:
    ndim, n, ite = 60, 47, 31
    xi  = vec([0.0000000000000000       0.44074620753680366        1.2190552270904040        2.3347760331235596        3.7887094567541970
        5.5718945978404708        7.6689972125944621        10.060434488482580        12.733770402821520        15.671866977274471
        18.842564947775848        22.211125483849631        25.740641931345998        29.392444169467371        33.126522288704130  
        36.901964947775845        40.677407606847567        44.411485726084329        48.063287964205706        51.592812898021521
        54.963348718681608        58.142735530561986        61.098971670733697        63.800414794198794        66.216309192726627
        68.317498170117346        70.077257917641631        71.472178009771781        72.483011008870633        73.095418495955613
        73.297402905362347        73.495381125374095        74.086800998417502        75.065298459040292        76.420152884300123
        78.136520205592305        80.195595543371837        82.574819237348777        85.248124014847377        88.186220589300319
        91.356918559801699        94.725479095875485        98.254995543371834        101.90679778149321        105.64087590072997
        109.41631855980170        113.19176121887341        0.0000000000000000        0.0000000000000000        0.0000000000000000  
        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000
        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000  ]') 
        
    bi  = vec([0.0000000000000000        2.4745511871254950        5.6535415917897565        9.0217956483937876        12.336388130013924        15.374699459348015        17.872596591715805        19.389265103084938        19.470334629888100        19.470334629888100
        19.470334629888100        19.470334629888100        19.470334629888100        19.470334629888100        19.470334629888100        19.470334629888100        19.470334629888100        19.470334629888100        19.470334629888100        19.421703717200348   
        18.696774608355046        17.220896249239004        15.166375014370439        12.722898324790991        10.088149927191735        7.4581743328295298        5.0179224470423360        2.9324013010399628        1.3388230888588288       0.34010363924408998  
        8.5025909811022496E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002 
        4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   4.2512954905511248E-002   0.0000000000000000        0.0000000000000000        0.0000000000000000  
        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000]')
        
    rni = vec([0.31533826832660750       0.57173742702119768       0.81219724788257108       0.90269258705773681       0.94672640854476298       0.97226522293211493       0.98894104468671684       0.99845196426994509       0.99999680694540605        1.0000000000000000
            1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000       0.99999937826068330       0.99982731522838886 
            0.99850586013772702       0.99567279498721784       0.99176744698370789       0.98718031684371466       0.98230291252748847       0.97750404726223072       0.97311231754692340       0.96940456031110123       0.96659922066293325       0.97648458361438661  
            0.99311284771759945       0.99967206764347949        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000   
            1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        1.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000   
            0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000]')  
            
    uinv =vec([0.0000000000000000       0.41232788144791949       0.69484769448346606       0.87101389669788887       0.98128982294812661        1.0572777365191697        1.1155410538547457        1.1478334518053532        1.0397120834699256        1.0208483066788829 
        1.0127227677468640        1.0089122361149454        1.0069272130443903        1.0059647578243049        1.0057023491673989        1.0060860938835354        1.0073117470299502        1.0101187892512191        1.0166872541232435        1.0522914806404915   
        1.0579707508297382        1.0458008481139760        1.0251962472333744       0.99887261065736377       0.96838145651404639       0.93415127095380235       0.89482945664006042       0.84503460603174330       0.76665868462116227       0.57832861017479698
        0.67956954870494446       0.91085945297732340       0.94987185259890972       0.96867306763204752       0.97940819877458840       0.98600310693451965       0.99023015937220771       0.99302181211083229       0.99491017491688016       0.99621419067431638  
        0.99713156304792716       0.99778802355998475       0.99826525037457403       0.99861730269630800       0.99888057062100188       0.99907993608762147       0.99923266171096192        0.0000000000000000        0.0000000000000000        0.0000000000000000 
        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000]')


    Reyn = 7362062.6848912220 
    Mach = 0.84 
    fexcr = 1.03  

    println("Benchmarking... _axisymm_BL")

    b = @benchmarkable aerodynamics._axisymm_BL($ndim, $n, $ite, $xi, $bi, $rni, $uinv,
    $Reyn, $Mach, $fexcr) seconds=30 evals=5
    bench_axisymm_BL = run(b)

    println("Benchmarking... fuselage_drag!")
    b = @benchmarkable aerodynamics.fuselage_drag!($(ac.fuselage), $parm, $para, $ipcruise1) seconds=30 evals=5
    bench_fuselage_drag = run(b)

    println("Benchmark results...")

    println("---------------------------------------")
    println("_axisymm_flow (FORTRAN on MacPro M2 ~ 30 μs)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"),bench_axisymm_flow)
    println(" ")

    println("---------------------------------------")
    println("_axisymm_BL (FORTRAN on MacPro M2 ~ 1.9 ms)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"),bench_axisymm_BL)
    println(" ")

    println("---------------------------------------")
    println("fuselage_drag! (FORTRAN on MacPro M2 ~ 1.95 ms)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"),bench_fuselage_drag)
    println(" ")
end

function benchmark_drag()
    # ====================
    #  Drag calculations
    # ====================

    println("\n---------------------------------------")
    println("Drag calculations")
    println("---------------------------------------")

    Re = 2e6
    println("Benchmarking... cfturb")
    bench_cfturb = @benchmark aerodynamics.cfturb($Re)


    println("Benchmarking... induced_drag!")
    bench_induced_drag = @benchmark aerodynamics.induced_drag!($view(para, :, ipcruise1), 
    $(ac.wing), $(ac.htail))

    nsurf = 2
    npout = [20, 10]
    npinn = [6, 0]
    npimg = [3, 2]
    Sref = 124.68530761144433
    bref =  35.486921631434697
    b = [35.486921631434697, 15.958117796995291]
    bs = [10.113772664958887, 1.5240000000000000]
    bo = [3.6067999999999998, 1.5240000000000000]
    bop = [0.72136000000000000, 1.5240000000000000]
    zcent = [-1.6764000000000001,  0.0000000000000000]
    po = [1.0000000000000000, 1.0000000000000000]
    gammat = [0.14999999999999999,  0.25000000000000000]
    gammas = [0.77000000000000002,  1.0000000000000000]
    fLo = -0.29999999999999999 
    ktip = 16
    specifies_CL = true
    CLsurfsp =[1.2502595056643222, 1.1976021933848557E-002] 
    idim::Int = 360
    jdim::Int = 360
    t     = zeros(Float64, jdim)
    y     = zeros(Float64, jdim)
    yp    = zeros(Float64, jdim)
    z     = zeros(Float64, jdim)
    zp    = zeros(Float64, jdim)
    gw    = zeros(Float64, jdim)

    yc    = zeros(Float64, idim)
    ycp   = zeros(Float64, idim)
    zc    = zeros(Float64, idim)
    zcp   = zeros(Float64, idim)
    gc    = zeros(Float64, idim)
    vc    = zeros(Float64, idim)
    wc    = zeros(Float64, idim)
    vnc   = zeros(Float64, idim)

    airfoil_section = aerodynamics.airtable(joinpath(TASOPT.__TASOPTroot__,"airfoil_data/C.air"));

    clp =  0.45     
    toc  = 0.126   
    Mperp  = 0.53  
    println("Benchmarking... airfun")
    bench = @benchmarkable aerodynamics.airfun($clp, $toc, $Mperp, 
                                            $airfoil_section) seconds = 30 evals = 100
    bench_airfun = run(bench)

    println("Benchmarking... _trefftz_analysis")
    bench = @benchmarkable aerodynamics._trefftz_analysis($nsurf, $npout,
                            $npinn, $npimg, 
                            $Sref, $bref,
                            $b,$bs,$bo,$bop, $zcent,
                            $po,$gammat,$gammas, $fLo, $ktip,
                        $specifies_CL,$CLsurfsp,
                            $t, $y, $yp, $z, $zp, $gw, $yc, $ycp, $zc, $zcp, $gc, $vc, $wc, $vnc) seconds=30 evals=100
    bench_trefftz = run(bench)

    wa = aerodynamics.tfp_workarrays(t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)
    tf = aerodynamics.tfp(nsurf, npout, npinn, npimg, 
    Sref, bref, b, bs, bo, bop, zcent, po, gammat, gammas, CLsurfsp, ktip, wa)
    bench = @benchmarkable aerodynamics.trefftz($tf) seconds=30 evals=100
    bench_trefftz_struct = run(bench)

    # You can use the following tempalte to profile code in the REPL
    # using ProfView
    # Profile.init(delay = 1e-6)
    # a = ProfileView.@profile (for i=1:10000; aerodynamics._trefftz_analysis(nsurf, npout,
    # npinn, npimg, 
    # Sref, bref,
    # b,bs,bo,bop, zcent,
    # po,gammat,gammas, fLo, ktip,
    # specifies_CL,CLsurfsp, t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc); end)
    # # Profile.print()
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    para1 = view(ac.para, :, ipcruise1, 1)
    gammat = wing.outboard.λ * para1[iarclt]
    gammas = wing.inboard.λ * para1[iarcls]
    Mach = para1[iaMach]
    fduo = para1[iafduo]
    fdus = para1[iafdus]
    fdut = para1[iafdut]
    CL   = para1[iaCL]
    CLhtail = para1[iaCLh]*htail.layout.S/wing.layout.S
    Reunit = para1[iaReunit]
    Reco = Reunit*wing.layout.root_chord
    rkSunsw = 0.5
    aRexp  = para1[iaaRexp]
    fexcdw = para1[iafexcdw]
    println("Benchmarking... wing_profiledrag_direct")
    bench = @benchmarkable aerodynamics.wing_profiledrag_direct($(ac.wing), $gammat, $gammas,
        $Mach, $CL, $CLhtail, $Reco,
        $aRexp, $rkSunsw, $fexcdw,
        $fduo, $fdus, $fdut) seconds=30 evals=5

    wing_profiledrag_direct = run(bench)

    println("Benchmarking... aircraft_drag!")
    bench = @benchmarkable aerodynamics.aircraft_drag!($parg, 
    $view(para, :, 10),
    $view(pare,:, 10), $(ac.wing), $(ac.htail), $(ac.vtail),
    $(1)) seconds=30 evals=1
    
    bench = @benchmarkable aerodynamics.aircraft_drag!(ac, 1, 10, true) seconds=30 evals=1
    bench_aircraft_drag = run(bench)


    println("---------------------------------------")
    println("cfturb (FORTRAN on MacPro M2 ~ --- μs)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"),bench_cfturb)
    println(" ")

    println("---------------------------------------")
    println("induced_drag! (FORTRAN on MacPro M2 ~ 4.5 μs)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"), bench_induced_drag)
    println(" ")
    
    println("---------------------------------------")
    println("airfun (FORTRAN on MacPro M2 ~ 420 ns)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"),bench_airfun)
    println(" ")

    println("---------------------------------------")
    println("trefftz (FORTRAN on MacPro M2 ~ 4.6 μs)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"),bench_trefftz)
    println(" ")

    println("---------------------------------------")
    println("trefftz struct (FORTRAN on MacPro M2 ~ 4.6 μs)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"),bench_trefftz_struct)
    println(" ")

    println("---------------------------------------")
    println("wing_profiledrag_direct (FORTRAN on MacPro M2 ~ 2.01 μs)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"),bench_wing_profiledrag_direct)
    println(" ")

    println("---------------------------------------")
    println("aircraft_drag! (FORTRAN on MacPro M2 ~ 6.7 μs)")
    println("---------------------------------------")
    show(stdout, MIME("text/plain"),bench_aircraft_drag)

end


function benchmark_gasfuns()
    t = [   175.00e0, 200.00e0, 225.00e0, 250.00e0, 275.00e0, 300.00e0,
        325.00e0, 350.00e0, 375.00e0, 400.00e0, 450.00e0, 500.00e0,
        550.00e0, 600.00e0, 650.00e0, 700.00e0, 750.00e0, 800.00e0,
        850.00e0, 900.00e0, 950.00e0, 1000.00e0, 1050.00e0, 1100.00e0,
        1150.00e0, 1200.00e0, 1250.00e0, 1300.00e0, 1350.00e0, 1400.00e0,
        1500.00e0, 1600.00e0, 1700.00e0, 1800.00e0, 1900.00e0, 2000.00e0,
        2100.00e0, 2200.00e0, 2300.00e0, 2400.00e0, 2500.00e0, 2600.00e0,
        2700.00e0, 2800.00e0, 2900.00e0, 3000.00e0, 3500.00e0, 4000.00e0,
        4500.00e0, 5000.00e0, 5500.00e0, 6000.00e0]
    tl = [ 5.16479e0, 5.29832e0, 5.41610e0, 5.52146e0, 5.61677e0, 5.70378e0,
        5.78383e0, 5.85793e0, 5.92693e0, 5.99146e0, 6.10925e0, 6.21461e0,
        6.30992e0, 6.39693e0, 6.47697e0, 6.55108e0, 6.62007e0, 6.68461e0,
        6.74524e0, 6.80239e0, 6.85646e0, 6.90776e0, 6.95655e0, 7.00307e0,
        7.04752e0, 7.09008e0, 7.13090e0, 7.17012e0, 7.20786e0, 7.24423e0,
        7.31322e0, 7.37776e0, 7.43838e0, 7.49554e0, 7.54961e0, 7.60090e0,
        7.64969e0, 7.69621e0, 7.74066e0, 7.78322e0, 7.82405e0, 7.86327e0,
        7.90101e0, 7.93737e0, 7.97247e0, 8.00637e0, 8.16052e0, 8.29405e0,
        8.41183e0, 8.51719e0, 8.61250e0, 8.69951e0]
    cp = [ 1039.00e0, 1039.00e0, 1039.00e0, 1039.10e0, 1039.30e0, 1039.60e0,
        1040.20e0, 1041.00e0, 1042.20e0, 1043.80e0, 1049.00e0, 1056.00e0,
        1065.00e0, 1075.00e0, 1086.00e0, 1098.00e0, 1110.00e0, 1122.00e0,
        1134.00e0, 1146.00e0, 1157.00e0, 1167.00e0, 1177.00e0, 1187.00e0,
        1196.00e0, 1204.00e0, 1212.00e0, 1219.00e0, 1226.00e0, 1232.00e0,
        1244.00e0, 1254.00e0, 1263.00e0, 1271.00e0, 1278.00e0, 1284.00e0,
        1290.00e0, 1295.00e0, 1300.00e0, 1304.00e0, 1307.00e0, 1311.00e0,
        1314.00e0, 1317.00e0, 1320.00e0, 1323.00e0, 1333.00e0, 1342.00e0,
        1349.00e0, 1355.00e0, 1362.00e0, 1369.00e0]
    cpt = [ 0.00051e0, -0.00051e0, 0.00154e0, 0.00638e0, 0.00894e0, 0.01784e0,
        0.02770e0, 0.03935e0, 0.05488e0, 0.07713e0, 0.12347e0, 0.16099e0,
        0.19255e0, 0.20879e0, 0.23230e0, 0.24203e0, 0.23959e0, 0.23962e0,
        0.24193e0, 0.23265e0, 0.20745e0, 0.19753e0, 0.20243e0, 0.19275e0,
        0.16659e0, 0.16090e0, 0.14980e0, 0.13990e0, 0.13058e0, 0.11776e0,
        0.11229e0, 0.09307e0, 0.08542e0, 0.07524e0, 0.06361e0, 0.06033e0,
        0.05506e0, 0.04944e0, 0.04717e0, 0.03186e0, 0.03539e0, 0.03659e0,
        0.02826e0, 0.03036e0, 0.03030e0, 0.02844e0, 0.01722e0, 0.01667e0,
        0.01211e0, 0.01288e0, 0.01437e0, 0.01363e0]
    h = [ -0.12781e6, -0.10184e6, -0.07586e6, -0.04989e6, -0.02391e6, 0.00208e6,
        0.02808e6, 0.05409e6, 0.08013e6, 0.10620e6, 0.15851e6, 0.21113e6,
        0.26415e6, 0.31765e6, 0.37167e6, 0.42626e6, 0.48146e6, 0.53726e6,
        0.59366e6, 0.65067e6, 0.70825e6, 0.76635e6, 0.82495e6, 0.88405e6,
        0.94363e6, 1.00363e6, 1.06403e6, 1.12481e6, 1.18594e6, 1.24739e6,
        1.37119e6, 1.49611e6, 1.62197e6, 1.74868e6, 1.87613e6, 2.00424e6,
        2.13294e6, 2.26220e6, 2.39195e6, 2.52216e6, 2.65271e6, 2.78361e6,
        2.91486e6, 3.04641e6, 3.17826e6, 3.31041e6, 3.97465e6, 4.64341e6,
        5.31625e6, 5.99224e6, 6.67146e6, 7.35422e6]
    s = [ -0.55312e3, -0.41438e3, -0.29201e3, -0.18253e3, -0.08349e3, 0.00695e3,
        0.09019e3, 0.16730e3, 0.23916e3, 0.30647e3, 0.42969e3, 0.54056e3,
        0.64162e3, 0.73471e3, 0.82118e3, 0.90210e3, 0.97826e3, 1.05028e3,
        1.11866e3, 1.18382e3, 1.24608e3, 1.30569e3, 1.36287e3, 1.41785e3,
        1.47082e3, 1.52189e3, 1.57121e3, 1.61888e3, 1.66502e3, 1.70971e3,
        1.79513e3, 1.87574e3, 1.95204e3, 2.02446e3, 2.09337e3, 2.15908e3,
        2.22187e3, 2.28200e3, 2.33968e3, 2.39510e3, 2.44839e3, 2.49973e3,
        2.54926e3, 2.59710e3, 2.64337e3, 2.68817e3, 2.89293e3, 3.07152e3,
        3.23001e3, 3.37245e3, 3.50192e3, 3.62073e3]

        bench = @benchmarkable gas_N2($t1, $t, $tl, $cp, $cpt, $h, $s) seconds=30 evals=50
end

function benchmark_gas()
    function f(n)
        for T in rand(200.0:500.0, n)
        _,_,cp,_ = TASOPT.engine.gas_N2(T)
        end
    end
    @benchmark f($100)
end

# benchmark_fuselage_drag()
# benchmark_drag()