//-----------------------------------------------------------------------------
//  Filename: printSolverDetails.fnc
// 
//  Description: 
//      Helper function to report Sovler Details
//-----------------------------------------------------------------------------

void printSolverDetails () {
  cout <<"\n"<< switchDes;
  cout << "\nDesign Independents:\n" << solver.independentNames;
  cout << "----------------------";
  cout << "\nDesign Dependents:\n" << solver.dependentNames;
  cout << "|=========================================================|" << endl;
  cout << endl;
}
