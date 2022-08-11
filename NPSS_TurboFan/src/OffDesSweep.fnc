
	real Thrusts [96];
	real j;
	real Wf_diff;
	
	for (j = 0; j<96; j++)
	{
		Thrusts[j] = 1.0 - 0.01*j;
	}
	
	cout << "ThrustLevels : " <<Thrusts <<endl;
	
	int i;
	real FAR;
	real PZ_phi;
	CASE = 1;
	for(i = 0; i<Thrusts.entries(); i++)
	{
		Eng.Pset
		{
			setOption("switchLimitSet","STD");
			setOption("switchParm","THRUST");
			parm_in = BaseThrust*Thrusts[i];
		}
		autoSolverSetup();
		//printSolverDetails();
		run();
		FAR = Eng.Core.BrnPri.Wfuel/Eng.Core.F035.W;
		PZ_phi = FAR/FAR_des*PZ_phi_des;
		
		Wf_diff = 9999;
		run();

		DesPt.update();
		CASE++;
		}

	cout<< CASE << " off-design case(s) run"<<endl;
