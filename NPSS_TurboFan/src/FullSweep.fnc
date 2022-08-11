// This "function" sweeps through thurst settings

PLAStream {filename = "../output/"+THERMPACKAGE+"FullSweep.out"; }

int details = 0;

if (!exists("Thrusts")){ real Thrusts [] ;}

Thrusts = {}; // Make sure we clear the array before redefining
Thrusts = linspace(1,0.045,50);

cout<<strFmt("SLS max thrust = %.3f lbf", SLS_thrust)<<endl;
CASE = 1;
  for(i = 0; i<Thrusts.entries(); i++)
	{
		Eng.Pset
		{
			setOption("switchLimitSet","STD");
			setOption("switchParm","THRUST");
			parm_in = SLS_thrust*Thrusts[i];
		}
				// solver.maxIterations=2000;
				// solver.maxJacobians = 500;
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
					CASE++;
	}
	OffDes.display();
	cout<< CASE << " off-design case(s) run"<<endl;
	
