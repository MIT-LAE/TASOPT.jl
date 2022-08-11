//-----------------------------------------------------------------------------
//  Filename: printSolverDetails.fnc
// 
//  Description: 
//      Helper function to report Sovler Details
//-----------------------------------------------------------------------------

void printSolverDetails () {
  // cout <<"\n"<< switchDes;
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