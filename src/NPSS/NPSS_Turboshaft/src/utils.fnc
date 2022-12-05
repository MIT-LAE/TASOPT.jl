

void whichDir(){
cout<<"\nCurrently in dir: \n\t"<<getcwd()<<endl;
}

void whichThermo(){
	cout<<strFmt("\nCurrently using %s thermo package", THERMPACKAGE)<<endl;
}

// Simple lines
void hashLine(){
	cout<<"#######################################"<<endl;
}

void hLine(){
	cout<<"_______________________________________"<<endl;
}

void dashLine(){
	cout<<"---------------------------------------"<<endl;
}


// Message between lines
void MsgBlock1(string msg){
	hashLine(); cout<<msg<<endl; hashLine();
}

void MsgBlock2(string msg){
	dashLine(); cout<<msg<<endl; dashLine();
}

void restart(string restartFileName){
	parseFile(restartFileName);
}