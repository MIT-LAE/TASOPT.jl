"""
PMSM sizes the electric machine for the motor and generator 
    of the turbo electric powertrain

Inputs:
 - 

Outputs:
 - 


Originally based of code by W. Enders
"""
function PMSM(P::Float64, N::Float64, rAsp, τAg, rSplit, p )

    # Setup/ assumptions
        # Stator
            ratSd = 1/50 #Slot depression to height ratio hsd/hsd
            kpf   = 0.5  #kpf = syrat = Packing factor
            z = 3        # phases
            Nshrt = 1    # Number of short pitched slots
            wSd   = 1e-6 # Slot depression width 
        # Rotor
            hRS = 0.002 #[m] Thickness of retaining sleeve for Vtip = 200 m/s
        #Shaft
            ratShft = 3/4 # Ratio of inner to outer dia for hollow Shaft
        
        # Material properties
            Br = 
    # Rotational speed
        ω = 2π*N    
    # ------------------
    # Machine geometry
    # ------------------
        Q = P/ω
        Qmax = 1.5*Q

        rRot = dRot/2.0

        hAg = ratAg*dRot + hRS # Air-gap Thickness

        hM  = ratM*hAg

        # Electric frequency
            rRoti = dRot/2 - (hM + hRS)
            f = p * N
            Ω = f * 2π
            Vtip = (rRoti + 0.5*hM)*ω

            #Catch if Vtip exceeds some values
        
        # Calculate Slot geometry
            NS  = ratSp * p  # Number of slots
            hS  = ratSM * hM # Slot height
            hSd = ratSd * hS # Slot depression height
            wST = 2π/Ns * (rRoti + hM + hAg + hSd + 0.5*hS) # Slot pitch

            wT = wST * ratW # Tooth width
            wS = wST - wT

            κS = wS/(rRoti + hAg - hRS + hSd + 0.5*hS)
            wSi = κS * (rRoti + hAg + hM + hSd)

            δ = wS/wST

            # Stator back Iron height
                hSBI = rRot/ratSplit - (rRot + hAg + hRS + hSd + hS)

            # Slots/pole/phases
                m = NS/(2*p*z)
                NSz = NS/z

                NSfpc = floor(m*z)
                NSc   = NSfpc - Nshrt
                cSp   = NSc/NS

            # Winding lengths
                l_ax  =  π * cSp * (rRoti+hAg+hM+hSd+0.5*hS) * δ/sqrt(1-δ^2)
                l_cir = 2π * cSp
                
                lET   = 2. * sqrt(l_ax^2 + (l_cir/2)^2)
            # Minimum magnet skew angle to prevent Cogging Torque
                κMs = 1/p * 180./π
            
        # Calc correction factors
            α = π * NSc/NSfpc
            kp = sin(α/2)

            γ  = 2π * p /NS
            kb = sin(m*γ/2)/(m*sin(γ/2))

            kw = kp*kb

            # Skew factor (Jagiela2013)
            κMs_rad = κMs * π/180.
            ks = sin(p*κMs_rad/2) / (p*κMs_rad/2)

            # Magentic gap factor
            rs = rRoti + hM + hAg
            r2 = rRoti + hM
            kg = ((rRoti^(p-1))/(rs^(2*p) - rRoti^(2*p))) * 
                 (  (p/(p+1))            * (r2^(p+1) - rRoti^(p+1)) + 
                    (p/(p-1)) * rs^(2*p) * (rRoti^(1-p) - r2^(1-p))  )
            
        # Calcualte back EMF
            Kc = 1/(1-(1/((wST/wS)*((5*hAg/wS)+1))));
            hAgEff = Kc*hAg        # Effective Air-gap
            Cphi = (p*κM)/180  # Flux concentration factor
            K1 = 0.95              # Leakage factor
            Kr = 1.05              # Reluctance factor
            murec = 1.05           # Recoil permeability

            PC  = hM/(hAgEff*Cphi);  # Permeance coefficient
            BAg = ((K1*Cphi)/(1+(Kr*murec/PC)))*Br;

        # Calculate magnetic flux and internal voltage
            κMrad = κM*(π/180);

            BC1 = (4/π)*BAg*kg*sin(p*κMrad/2);
            λ = 2*rs*lRot*NSz*kw*ks*BC1/p;

            ErmsD  = omega*λ/sqrt(2); # RMS back voltage
    
        # Calculation of inductance and total reactance
            # Air-gap inductance
                LAg = (z)*(2/π)*(mu0*NSz^2*kw^2*lRot*rs)/(p^2*(hAg+hM));

            # Slot leakage inductance
                Cperm = mu0*((1/3)*(hS/wSi) + hSd/wSi);
                Las = 2*p*lRot*Cperm*(4*(m-NShrt)+2*NShrt);
                Lam = 2*p*lRot*NShrt*Cperm;
                Ls = Las - Lam; # equ. for 3 phases only

            # End-turn inductance 
                areaS = wS*hS; # Slot area
                Le = 0.25*(mu0*wST*NSz^2)*log(wST*sqrt(π)/sqrt(2*areaS));
                #wST only true, if coils are placed in neighboring slots 


            # Total inductance and reactance per phase
                Ltot = LAg+Ls+Le;
                Xtot = omega*Ltot;


        ## ------------Calculation of machine dimensions and weights---------------
        #Armature
            #Total Armature length per phase (assuming two coils per slot)
                lArm=2*NSz*(lRot+lEt); 
            # Armature conductor area
                areaArm = 0.5*areaS*kpf;
            # Total mass of armature conductor #mass=pha*l*area*rho
                mArm = z*lArm*areaArm*rhoCon;
        
        
        #Iron /Stator Core
            # SBI inside radius
                rSBIi = rRoti+hM+hAg+hSd+hS;
            # SBI outside radius
                rSBIo = rSBIi+hSBI;
            # Core mass
                mSBI = pi*(rSBIo^2-rSBIi^2)*lRot*rhoIron; # SBI
                mTeeth = (NS*wT*hS+2*pi*(rRoti+hAg)*hSd-NS*hSd*wSd)*lRot*rhoIron; # Teeth
                mIron = mSBI + mTeeth;
        
        # Magnet mass
            mM = (p*kappaMrad)*((rRoti+hM)^2-rRoti^2)*lRot*rhoMag; 
        
        # Shaft mass (Hollow shaft)
            #tauMax=Qmax/Wt Wt=section modulus (Widerstandsmoment)
            #Wt=pi/16*(da^4-di^4)/da, (Dankert,Technische Mechanik, 2018)
            lShft  = 1.2*lRot; #Assumption to account for bearings, endcaps, etc.
            rShfto = 0.5*((16*Qmax)/(tauMax*pi*(1-ratShft^4)))^(1/3); #Outer shaft radius
            mShft  = (pi*rShfto^2*(1-ratShft^2))*lShft*rhoShft; #Mass of shaft
            tShft  = rShfto*(1-ratShft); #thickness of shaft, just for comparison
        
        # Service mass fraction
            #kServ = 1.15; # Rucker2005
            #kServ= 1.5;   # Yoon2016 (Tab.4) Excluding Ground Cylinder and Heat Sink
            kServ= 1.7;    # Yoon2016 (Tab.4) Including Ground Cylinder and Heat Sink
        
        
        # Total mass
        mPMSM = kServ*(mIron+mShft+mM+mArm); 
        
        
        # Final machine dimensions
        lPMSM = lRot+2*lEtAx; #Total length
        dPMSM = 2*rSBIo; # Total diameter without housing



end
