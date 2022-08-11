//Create intermediary points given 2 arrays

real[] linspace(real start, real end, int n) {
	int i;
	real step = (end-start)/(n-1);
	real arr[]={};
	for(i= 0; i<n; i++)
	{
		arr.append(start+i*step);
	}
	return arr;
}

//Test
// cout<<"testing..."<<endl;
// real test[] = linspace(1,10, 10);
// cout<<test<<endl;