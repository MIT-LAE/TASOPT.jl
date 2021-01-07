void SaveIndepVals (string Sname, string out){
	
	 if (exists("saveStream")){
		saveStream.filename = out; 
	 }
	 else{
		OutFileStream saveStream{filename = out; }
	 }
	// cout<<"Saving model state to : "<<out<<endl;
	saveStream <<"//"<< Sname << " Independent Variables:" << endl;
	saveStream<<"solver.independentNames = "<<Sname->independentNames<<";"<<endl;
	saveStream<<"solver.independentValues = "<<Sname->independentValues<<";"<<endl;
	
}