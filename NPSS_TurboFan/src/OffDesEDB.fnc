
PLAStream {filename = "../output/OffDesEDB.out"; }

	if(EngineName == "CFM56-5B"){
	//Thrust ratios from EDB for the CFM 56_5B---------------------------------------------------
		real SLS_thrust = 32000;
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
		 real SLS_thrust = 27293;
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
	real FAR100;
	real FAR7;
	real PZ_phi;
	real RMS_error=0;
	
	real LPC_PR100, LPC_PR7, LPC_eff100, LPC_eff7, LPC_R100, LPC_R7, LPC_Nc100, LPC_Nc7;
	real FAR100q7;
	real phi100, phi7;
	
	phi100 = 1.85;
	
	Eng.Amb.alt_in = 0;
	Eng.Amb.MN_in  = 0;
	
	CASE = 1;
	for(i = 0; i<Thrusts.entries(); i++)
	{
		Eng.Pset
		{
			setOption("switchLimitSet","STD");
			setOption("switchParm","THRUST");
			parm_in = SLS_thrust*Thrusts[i];
		}
		Eng.Core.CmpH.switchAud = "AUDIT";
		Eng.Core.CmpH.s_effAud = 1;
		if (i>18){
			Eng.Core.CmpH.s_effAud = 1.0;
		}
		autoSolverSetup();
		//printSolverDetails();
		run();
		if(solver.converged == 0)
		{
			solver.defaultDxLimit = 0.05;
			autoSolverSetup();
			run();
		}
		FAR = Eng.Core.BrnPri.FAR;
		if(i==0) {
			FAR100 = FAR;
			LPC_PR100 = Eng.CmpL.PR;
			LPC_eff100 = Eng.CmpL.eff;
			LPC_R100   = Eng.CmpL.S_map.RlineMap;
			LPC_Nc100  = Eng.CmpL.S_map.NcMap;
		}
		PZ_phi = FAR/FAR_des*PZ_phi_des;
		
		Wf_diff = 100*(Wf_EDB[i] - convertUnits("Eng.Core.BrnPri.Wfuel","kg/sec"))/Wf_EDB[i];
		//cout<<Wf_diff<<endl;
		RMS_error += (1 - convertUnits("Eng.Core.BrnPri.Wfuel","kg/sec")/Wf_EDB[i] )**2;
		//cout<<RMS_error<<endl;
		//run();

		OffDes.update();
		CASE++;
	}

	FAR7 = FAR;
	FAR100q7 = FAR100/FAR7;
	LPC_PR7 = Eng.CmpL.PR;
	LPC_eff7 = Eng.CmpL.eff;
	LPC_R7   = Eng.CmpL.S_map.RlineMap;
	LPC_Nc7  = Eng.CmpL.S_map.NcMap;
	
	OffDes.display();
	cout<< CASE << " off-design case(s) run"<<endl;
	
	RMS_error = 100*(sqrt(RMS_error/Wf_EDB.entries()));
	cout<<"------RMS error for all points = "<<RMS_error<<" %------\n"<<endl;
	cout<<"\tFAR100/FAR7 = "<<FAR100q7<<endl<<endl;
	cout<<"\tFAR100 = "<< FAR100<<"\n\tFAR7 = "<<FAR7<<endl;
	
	cout<<"\n\tLPC    100%Foo\t7%Foo"<<endl;
	cout<<"\t---------------------"<<endl;
	cout<<strFmt("\tPR     %.3f\t%.3f", LPC_PR100, LPC_PR7)<<endl;
	cout<<strFmt("\teff    %.3f\t%.3f",LPC_eff100, LPC_eff7)<<endl;
	cout<<strFmt("\tR      %.3f\t%.3f" ,LPC_R100, LPC_R7)<<endl;
	cout<<strFmt("\tNc     %.3f\t%.3f" ,LPC_Nc100, LPC_Nc7)<<endl;
	cout<<strFmt("\tphi    %.3f\t%.3f" ,phi100, phi100/FAR100q7)<<endl<<endl;
