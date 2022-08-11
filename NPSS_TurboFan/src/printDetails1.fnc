//-----------------------------------------------------------------------------
//  Filename: printSolverDetails.fnc
// 
//  Description: 
//      Helper function to report Sovler Details
//-----------------------------------------------------------------------------

void printRunDetails () {
	
	if(Eng.switchDes == "DESIGN")
	{
	cout << "   Design Case " << CASE << ", ESIseverityMax = " << ESIseverityMax << endl;

	cout << "\n---------------- SUMMARY --------------------" << endl;	
	
	cout << "   Alt = "<<Eng.Amb.alt_in << "      MN = "<<Eng.Amb.MN_in <<endl;
	cout << " OPR " << Eng.Cycle.OPR << ", Wf = "<< Eng.Core.BrnPri.Wfuel<<" "+toStr(Eng.Core.BrnPri.Wfuel.units)<<", SFC = "<< Eng.Perf.SFC<< endl;

	cout <<  "D E S I G N   C A S E - " << CASE << endl;
	cout <<  "BPR = "<< Eng.SpltFan.BPR << endl;
	cout <<  "FPR = "<< Eng.CmpF.PR << endl;
	cout <<  "LPC PR = "<< Eng.CmpL.PR << endl;
	cout <<  "HPC PR = "<< Eng.Core.CmpH.PR << endl;
	cout <<  "OPR = "<< Eng.Cycle.OPR << endl;
	cout <<  "\nTt40[R] = "<<Eng.Core.BrnPri.TtCombOut<< endl;
	cout <<  "Tt41[R] = " << Eng.Core.TrbH.F041.Tt << endl; // REMEMBER THAT the turbine has a new port F041 that represents cooled flow!
	cout <<  "Fan Dia [in] = "<<  (Eng.CmpF.Fl_I.Aphy*4/3.1415/(1-0.34**2))**0.5 << endl;
	cout <<  "TrbH_eff = " << Eng.Core.TrbH.effDes << endl; 
	cout <<  "TrbL_eff = " << Eng.TrbL.effDes << endl; 
	
	cout <<  "---------------- M A P   P O I N T S --------------------" << endl;
	cout <<  "CmpF [s_PR,sWp,Seff]:\t" << Eng.CmpF.S_map.s_PRdes<<", "<<Eng.CmpF.S_map.s_WcDes<<", "<< Eng.CmpF.S_map.s_effDes << endl;
	cout <<  "CmpL [s_PR,sWp,Seff]:\t" << Eng.CmpL.S_map.s_PRdes<<", "<<Eng.CmpL.S_map.s_WcDes<<", "<< Eng.CmpL.S_map.s_effDes << endl;
	cout <<  "CmpH Map Scalars [s_PR,sWp,Seff,sNc]:\t" << Eng.Core.CmpH.S_map.s_PRdes<<", "
	<<Eng.Core.CmpH.S_map.s_WcDes<<", "<< Eng.Core.CmpH.S_map.s_effDes <<", "<< Eng.Core.CmpH.S_map.s_NcDes << endl;
	
	cout <<  "CmpH Map POINTS [PR, Wc,eff,Nc,RlineMap]:\t" << Eng.Core.CmpH.S_map.PRmap<<", "
	<<Eng.Core.CmpH.S_map.WcMap<<", "<< Eng.Core.CmpH.S_map.effMap <<", "<< Eng.Core.CmpH.S_map.NcMap <<", "<< Eng.Core.CmpH.S_map.RlineMap << endl;

	cout <<  "\nTrbH [s_PR,sWp,Seff,sNp]:\t" << Eng.Core.TrbH.S_map.s_PRdes<<", "
	<<Eng.Core.TrbH.S_map.s_WpDes<<", "<< Eng.Core.TrbH.S_map.s_effDes <<", "<< Eng.Core.TrbH.S_map.s_NpDes << endl;
	cout <<  "TrbL [s_PR,sWp,Seff]:\t" << Eng.TrbL.S_map.s_PRdes<<", "<<Eng.TrbL.S_map.s_WpDes<<", "<< Eng.TrbL.S_map.s_effDes << endl;
	cout <<  "-----------------------------------------------------------"<<endl;
	cout << endl;

	cout << "--St F020: Fan face--"<<endl;
	cout << "Tt2 = "<<Eng.F020.Tt << " R = "<< convertUnits("Eng.F020.Tt","K")<< " K"<< endl;
	cout << "Pt2 = "<<Eng.F020.Pt<< " psi = "<< convertUnits("Eng.F020.Pt","kPa")<< " kPa"<< endl;
	cout << "W2 = "<<Eng.F020.W<< " lbm/s = "<< convertUnits("Eng.F020.W","kg/sec")<< " kg/s"<< endl;
	cout << "\n--St F035: Combustor inlet--"<<endl;
	cout << "Tt3 = "<<Eng.Core.F035.Tt << " R = "<< convertUnits("Eng.Core.F035.Tt","K")<< " K"<< endl;
	cout << "Pt3 = "<<Eng.Core.F035.Pt<< " psi = "<< convertUnits("Eng.Core.F035.Pt","kPa")<< " kPa"<< endl;
	cout << "W3 = "<<Eng.Core.F035.W<< " lbm/s = "<< convertUnits("Eng.Core.F035.W","kg/sec")<< " kg/s"<< endl;
	cout << "\n--St F040: Combustor exit--"<<endl;
	cout << "Tt40 = "<<Eng.Core.F040.Tt << " R = "<< convertUnits("Eng.Core.F040.Tt","K")<< " K"<< endl;
	cout << "Pt40 = "<<Eng.Core.F040.Pt<< " psi = "<< convertUnits("Eng.Core.F040.Pt","kPa")<< " kPa"<< endl;
	cout << "W40 = "<<Eng.Core.F040.W<< " lbm/s = "<< convertUnits("Eng.Core.F040.W","kg/sec")<< " kg/s"<< endl;
	cout <<  "-----------------------------------------------------------"<<endl;

	
	cout <<  "Fuel flow Rate = "<< Eng.Core.BrnPri.Wfuel<< " lbm/s = "<< convertUnits("Eng.Core.BrnPri.Wfuel","kg/sec")<< " kg/s"<< endl;
	cout <<  "Net Thrust [lbf] = "<< Eng.Perf.Fn << " lbf = "<< convertUnits("Eng.Perf.Fn","kN")<< " kN"<< endl;
	cout <<  "SFC [lbm/hr/lbf] = "<< Eng.Perf.SFC << endl;
	cout <<  "Core Design status: "<<Eng.Core.TrbH.switchDes << endl;
	cout << " ------ COOLING FLOWS -------"<<endl;
	cout << "Non Chargeable Flow   : "<<Eng.Core.TrbH.TCLA_NC.W<<" "<<Eng.Core.TrbH.TCLA_NC.W.units
	<< " ---" <<100*Eng.Core.TrbH.TCLA_NC.W/Eng.Core.F030.W<< "%"<<endl;
	cout << "Chargeable Flow       : "<<Eng.Core.TrbH.TCLA_CH.W<<" "<<Eng.Core.TrbH.TCLA_CH.W.units
	<< " ---" <<100*Eng.Core.TrbH.TCLA_CH.W/Eng.Core.F030.W<< "%"<<endl;
	cout << "Fraction of core flow : "<<1-Eng.Core.F035.W/Eng.Core.F030.W<<endl;
	
	cout<< CASE << " design case(s) run"<<endl;
	cout<<"--------------------------------------------"<<endl;
	endl;
	}
	else
	{
	cout << "----------------------------"<<endl;
	cout << "T H R U S T - L E V E L  =  "<<Thrusts[i]<<endl;
	cout <<  "BPR = "<< Eng.SpltFan.BPR << endl;
	cout <<  "FPR = "<< Eng.CmpF.PR << endl;
	cout <<  "OPR = "<< Eng.Cycle.OPR << endl;
	cout <<  "Tt4 = " << Eng.Core.BrnPri.Fl_O.Tt << " R = "<< convertUnits("Eng.Core.BrnPri.Fl_O.Tt","K")<< " K"<< endl;
	
	cout << "Tt3 = "<<Eng.Core.F035.Tt << " R = "<< convertUnits("Eng.Core.F035.Tt","K")<< " K"<< endl;
	cout << "Pt3 = "<<Eng.Core.F035.Pt<< " psi = "<< convertUnits("Eng.Core.F035.Pt","kPa")<< " kPa"<< endl;
	cout << "W3 = "<<Eng.Core.F035.W<< " lbm/s = "<< convertUnits("Eng.Core.F035.W","kg/sec")<< " kg/s"<< endl;
	cout <<  "Fuel flow Rate = "<< Eng.Core.BrnPri.Wfuel<< " lbm/s = "<< convertUnits("Eng.Core.BrnPri.Wfuel","kg/sec")<< " kg/s"<< endl;
	cout <<  "Net Thrust [lbf] = "<< Eng.Perf.Fn << " lbf = "<< convertUnits("Eng.Perf.Fn","kN")<< " kN"<< endl;
	cout << "----------------------------"<<endl;
	}
	

}

void printSolverDetails () {
  cout <<"\n"<< switchDes;
  cout << "\nDesign Independents:\n" << solver.independentNames;
  cout << "----------------------";
  cout << "\nDesign Dependents:\n" << solver.dependentNames;
  cout << "|=========================================================|" << endl;
  cout << endl;
}

void printSolverConditions (){
		string Sname = "solver";
		string indeps[] = Sname->independentNames;
		real indepVals[] = Sname -> independentValues;

		string deps[] = Sname->list("Dependent", FALSE );
		cout << Sname << " Independent Variables:" << endl;
		int i;
		string var;
		for (i=0; i<indeps.entries(); i++) {
		var = indeps[i]->varName;
		cout << " " << i+1 << " " << indeps[i] << ": \n\t"
		<< var <<" = "<<indepVals[i]<< var.units << endl;
		}
		
		cout << Sname << " Dependent Conditions:" << endl;
		string lhs, rhs;
		for (i=0; i<deps.entries(); i++) {
		lhs = deps[i]->eq_lhs;
		rhs = deps[i]->eq_rhs;
		cout << " " << i+1 << " " << deps[i] << ": \n\t"
		<< lhs << " = " << rhs << endl;
		}

}