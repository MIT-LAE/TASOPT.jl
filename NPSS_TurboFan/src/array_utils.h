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

real[] logspace(real start, real end, int n) {
	int i;
	real base = 10;
	real step = (end-start)/(n-1);
	real arr[]={};
	for(i= 0; i<n; i++)
	{
		arr.append(base**(start+i*step));
	}
	return arr;
}

// //Test
// cout<<"testing..."<<endl;
// real test[] = linspace(1,10, 10);
// cout<<test<<endl;
// cout<<"testing..."<<endl;
// real test2[] = logspace(2,3, 4);
// cout<<test2<<endl;