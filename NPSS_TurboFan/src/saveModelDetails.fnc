

// Output map scalars:
OutFileStream outfile {filename = "../mapScalars_7B.out"; }
void saveMapScalars() {
	outfile<<"//[prash] High pressure compressor scalars calculated by CFM56-7B27 model"<<endl;
	outfile<<"Eng.Core.CmpH.S_map.s_PRdes = "<<Eng.Core.CmpH.S_map.s_PRdes<<";"<<endl;
	outfile<<"Eng.Core.CmpH.S_map.s_WcDes = "<<Eng.Core.CmpH.S_map.s_WcDes<<";"<<endl;
	outfile<<"Eng.Core.CmpH.S_map.s_NcDes = "<<Eng.Core.CmpH.S_map.s_NcDes<<";"<<endl;
	outfile<<"Eng.Core.CmpH.S_map.s_effDes = "<<Eng.Core.CmpH.S_map.s_effDes<<";"<<endl<<endl;
	
	outfile<<"//[prash] High pressure Turbine scalars calculated by CFM56-7B27 model"<<endl;
	outfile<<"Eng.Core.TrbH.S_map.s_PRdes = "<<Eng.Core.TrbH.S_map.s_PRdes<<";"<<endl;
	outfile<<"Eng.Core.TrbH.S_map.s_WpDes = "<<Eng.Core.TrbH.S_map.s_WpDes<<";"<<endl;
	outfile<<"Eng.Core.TrbH.S_map.s_effDes = "<<Eng.Core.TrbH.S_map.s_effDes<<";"<<endl;
	outfile<<"Eng.Core.TrbH.S_map.s_NpDes = "<<Eng.Core.TrbH.S_map.s_NpDes<<";"<<endl<<endl;
	
	outfile<<"//[prash] Bleed flow ratios calculated by CFM56-7B27 model"<<endl;
	outfile<<"Eng.Core.TrbH.TCLA_NC.fracW = "<<Eng.Core.TrbH.TCLA_NC.fracW<<";"<<endl;
	outfile<<"Eng.Core.TrbH.TCLA_CH.fracW = "<<Eng.Core.TrbH.TCLA_CH.fracW<<";"<<endl<<endl;
	
	outfile<<"//[prash] Nc design calculated by CFM56-7B27 model"<<endl;
	outfile<<"Eng.Core.CmpH.NcDes = "<<Eng.Core.CmpH.NcDes<<";"<<endl;
	outfile<<"Eng.Core.TrbH.NpDes = "<<Eng.Core.TrbH.NpDes<<";"<<endl;
	
}

// Output polytropic efficiencies:
OutFileStream polytropic{filename = "../effPolytropic_7B.out"; }
void saveEffPoly(){
	polytropic<<"real Fan_effPoly = "<<Eng.CmpF.effPoly<<";"<<endl;
	polytropic<<"real LPC_effPoly = "<<Eng.CmpL.effPoly<<";"<<endl;
	polytropic<<"real HPC_effPoly = "<<Eng.Core.CmpH.effPoly<<";"<<endl;
	polytropic<<"real HPT_effPoly = "<<Eng.Core.TrbH.effPoly<<";"<<endl;
	polytropic<<"real LPT_effPoly = "<<Eng.TrbL.effPoly<<";"<<endl;
	
}

// [TODO] Output flow_areas
OutFileStream flow_areas{filename = "../flow_areas_7B.out"; }
void saveAreas(){
	flow_areas<<"Core.CmpH.Fl_O.Aphy = "<<Eng.Core.CmpH.Fl_O.Aphy<<";"<<endl;
}

void saveModelDetails(){
	saveMapScalars();
	saveEffPoly();
}