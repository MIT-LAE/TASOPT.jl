
// Converts a MAP to python np arrays
void convertCompressorMap(string name, string MAPNAME){
	
	parseFile(MAPNAME);
	
	
	OutFileStream file {filename = name; }
	
	file {filename = name; }
	real ALPHA[] = S_map.TB_Wc.getValues("ALPHA");
	real SPED[]  = S_map.TB_Wc.getValues("SPED",ALPHA[0]);
	real R[] = S_map.TB_Wc.getValues("R",ALPHA[0],SPED[0]);
	cout<<"Converting "<< MAPNAME <<" to a numpy array ..."<<endl;
	cout<<"Alpha = "<<ALPHA<<"\nSpeed ="<<SPED<<"\nRline ="<<R<<endl;
	int i;
	int j;
	
	
	file<<"# MAP NAME : "<<MAPNAME<<endl<<endl;
	file<<"# ---------------------------------"<<endl;
	file<<"import numpy as np"<<endl<<endl;
	file<<"class CmpMap(object):"<<endl;
	file<<"    pass"<<endl;
	file<<"CmpMap.NcDes = "<<S_map.NcMapDes<<endl;
	file<<"CmpMap.RlineDes = "<<S_map.RlineMapDes<<endl;

	file<<"\nCmpMap.Nc = np.array([";
	for(i =0; i<SPED.entries(); i++){ file<<SPED[i]<<",";}
	file<<"])\n\n";
	file<<"\nCmpMap.R = np.array([";
	for(i =0; i<R.entries(); i++){ file<<R[i]<<",";}
	file<<"])\n\n";



	file<<"\nCmpMap.Wc = np.array([";
	for(i =0; i<SPED.entries(); i++){
		file<<"[";
		for(j=0; j<R.entries(); j++){
			
			file<<S_map.TB_Wc(0,SPED[i],R[j])<<",";
			
		}
		file<<"],\n";
	}
	file<<"])"<<endl<<endl;

	file<<"\nCmpMap.eff = np.array([";
	for(i =0; i<SPED.entries(); i++){
		file<<"[";
		for(j=0; j<R.entries(); j++){
			
			file<<S_map.TB_eff(0,SPED[i],R[j])<<",";
			
		}
		file<<"],\n";
	}
	file<<"])"<<endl<<endl;

	file<<"\nCmpMap.PR = np.array([";
	for(i =0; i<SPED.entries(); i++){
		file<<"[";
		for(j=0; j<R.entries(); j++){
			
			file<<S_map.TB_PR(0,SPED[i],R[j])<<",";
			
		}
		file<<"],\n";
	}
	file<<"])"<<endl<<endl;
	file.close();
	delete("S_map"); // use this to delete objects!

}

// convertCompressorMap("lpcCFM56.py", "lpcCFM56.map");
convertCompressorMap("fanCFM56.py", "fanCFM56.map");
// convertCompressorMap("hpcCFM56.py", "hpcCFM56.map");
// convertCompressorMap("lpcE3.py", "lpcE3.map");
// convertCompressorMap("hpcE3.py", "hpcE3.map");


// convertCompressorMap("CFM56_padded.py", "lpcCFM56_prash.map");
convertCompressorMap("CFM56_padded.py", "lpcCFM56_prash_new.map");

// cout<<S_map.TB_Wc.getTextRep(0)<<endl;
// cout<<S_map.TB_Wc.getIndependentNames()<<endl;


