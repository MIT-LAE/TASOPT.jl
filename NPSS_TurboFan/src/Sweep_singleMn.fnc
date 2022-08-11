// This "function" sweeps through differnet PLA's 
// Note that the loop zig-zags to help with convergnece

  
int details  = 0;

real PLAs [] = {120,110,100,98,95,93,90,88,85,83,80,75,70,65,60,55,50,40,30,20,10};

int plaLoop; int altLoop;
int i=-1;


//Eng.Amb.alt_in=40000;

cout<<"Alt = " << Eng.Amb.alt_in<<endl;

  for(plaLoop=0;plaLoop<PLAs.entries();plaLoop++)
	{

			Eng.Pset
			{
				setOption("switchLimitSet","STD");
				setOption("switchParm","PLA");
				parm_in =PLAs[plaLoop];
			}
			solver.maxIterations=1000;
			solver.maxJacobians = 200;
			autoSolverSetup();
			solver.forceNewJacobian=TRUE;
			
			if (details == 1) { printSolverDetails(); }
			run();
			SSTF_out.update();
			cout<<"CASE: "<< CASE << "; Mn0 = "<<Eng.Amb.MN_in<< endl;
			CASE++;

	} 
	

