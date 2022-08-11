// This "function" sweeps through differnet PLA's 
// Note that the loop zig-zags to help with convergnece

real Mn_max  = 2.11;
real Mn_min  = 0.09;
real Mn_step = 0.4;
  
int details  = 0;

real Thrusts [] = {5000,5500,6000,6500,7000,7500,8000,9000,10000,15000};
real Alts [] = {40000.0, 20000.0, 0.0};

int thrustLoop; int altLoop;
int i=-1;


//Eng.Amb.alt_in=40000;

cout<<"Alt = " << Eng.Amb.alt_in<<endl;

  for(thrustLoop=0;thrustLoop<Thrusts.entries();thrustLoop++)
	{

			Eng.Pset
			{
				setOption("switchLimitSet","NONE");
				setOption("switchParm","THRUST");
				parm_in =Thrusts[thrustLoop];
			}
			solver.maxIterations=100;
			solver.maxJacobians = 200;
			autoSolverSetup();
			solver.forceNewJacobian=TRUE;
			
			if (details == 1) { printSolverDetails(); }
			run();
			SSTF_out.update();
			
			cout<<"CASE: "<< CASE << "; Mn0 = "<<Eng.Amb.MN_in<< endl;
			CASE++;
	
	} 
	

