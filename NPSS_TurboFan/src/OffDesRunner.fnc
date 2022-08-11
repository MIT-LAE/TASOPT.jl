void OffDes_EDB(string EngineName)
{
	cout<< EngineName<<endl;
	if(EngineName == "CFM56-5B"){
	//Thrust ratios from EDB for the CFM 56_5B---------------------------------------------------

		real Thrusts [] = { 1.  ,0.96907941, 0.93745608, 0.85, 0.84399157 ,0.8237175,
		 0.79683767, 0.73436402, 0.72803935, 0.71739283, 0.68798313, 0.6753338,
		 0.62420942, 0.61883345, 0.58478566, 0.57403373, 0.3,        0.29072382,
		 0.28123682, 0.25319747, 0.22030921, 0.21841181, 0.20639494, 0.20260014,
		 0.07,       0.06783556, 0.06562193, 0.05907941, 0.05096275, 0.04815882, 0.04727337};
		 
		real Wf_EDB [] = { 1.462, 1.385, 1.318, 1.153, 1.142, 1.107, 1.063, 0.965, 0.956, 0.939,
		 0.894, 0.875, 0.8  , 0.793, 0.743, 0.727, 0.369, 0.358, 0.347, 0.316,
		 0.279, 0.278, 0.264, 0.26 , 0.113, 0.111, 0.109, 0.102, 0.095, 0.092,
		 0.091 };  
	}
	else if(EngineName == "CFM56-7B"){
		 real Thrusts [] = { 1.        , 0.96375618, 0.88632619, 0.85      , 0.83196046, 0.81919275,
		 0.75453048, 0.75337727, 0.71416804, 0.70716639, 0.64135091, 0.60704283,
		 0.3       , 0.28912685, 0.26589786, 0.24958814, 0.22635914, 0.21425041,
		 0.07      , 0.06746293, 0.06204283, 0.05823723, 0.05281713, 0.04999176 };
		 
		 real Wf_EDB []= { 1.293, 1.213, 1.086, 1.031, 1.004, 0.986, 0.896, 0.895, 0.842, 0.832,
		 0.746, 0.702, 0.343, 0.331, 0.308, 0.291, 0.268, 0.256, 0.11 , 0.108,
		 0.103, 0.099, 0.094, 0.092 };
		
	}
	
	real Wf_diff;
	// --------------------------------------------------------- end EDB comparison
		   
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
		
		Wf_diff = 100*(Wf_EDB[i] - Eng.Core.BrnPri.Wfuel/2.205)/Wf_EDB[i];
		run();

		DesPt.update();
		CASE++;
		}

	cout<< CASE << " off-design case(s) run"<<endl;
}

void OffDes_sweep()
{
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
}