// This "function" sweeps through differnet PLA's 
// Note that the loop zig-zags to help with convergnece

real Mn_max  = 0.85;
real Mn_min  = 0.40;
real Mn_step = 0.1;

int details  = 0;

real PLAs [] = {110,106,104,102,100,98,95,93,90,88,85,83,80,75,70,65,60,55,50,40,30,20,10};
//real Alts [] = {40000.0, 20000.0, 0.0};

int plaLoop; int altLoop;
int i=-1;


//Eng.Amb.alt_in=40000;
//cout<<"Alt = " << Eng.Amb.alt_in<<endl;

  for(plaLoop=0;plaLoop<PLAs.entries();plaLoop++)
	{
		i=-1;
		i=i**(plaLoop+1);

		
		while(Eng.Amb.MN_in<=Mn_max && Eng.Amb.MN_in>=Mn_min )
		{
			Eng.Pset
			{
				setOption("switchLimitSet","STD");
				setOption("switchParm","PLA");
				parm_in =PLAs[plaLoop];
			}
			solver.maxIterations=2000;
			solver.maxJacobians = 500;
			autoSolverSetup();
			solver.forceNewJacobian=TRUE;
			
			if (details == 1) { printSolverDetails(); }
			run();
			if(solver.converged == 0)
			{
				solver.defaultDxLimit = 0.05;
				autoSolverSetup();
				run();
			}
			OffDes.update();
			cout<<"CASE: "<< CASE << "; Mn0 = "<<Eng.Amb.MN_in<< endl;
			CASE++;
			Eng.Amb.MN_in = Eng.Amb.MN_in + Mn_step*i;
	

		}
		Eng.Amb.MN_in = Eng.Amb.MN_in - Mn_step*i;
	} 
	

