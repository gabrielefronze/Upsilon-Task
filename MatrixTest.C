using namespace std;

void MatrixTest(){

	const int fullmatDim=13;
	const int dim1S=7;
	const int dim2S=8;
	const int dim3S=8;
	const int shift=4;

	int fullmat[fullmatDim*fullmatDim];
	for(int i=0;i<fullmatDim*fullmatDim; i++){
		fullmat[i]=i;
	}

	int Upsilon1Smat[dim1S*dim1S];
	for(int i=0;i<dim1S;i++){
		for(int j=0;j<dim1S;j++){
			Upsilon1Smat[dim1S*i+j]=fullmat[fullmatDim*shift+shift+j+fullmatDim*i];
		}
	}

	int Upsilon2Smat[dim2S*dim2S];
	for(int i=0;i<dim2S;i++){
		for(int j=0;j<dim2S;j++){
			Upsilon2Smat[dim2S*i+j]=fullmat[fullmatDim*shift+shift+j+fullmatDim*i];
		}
	}

	int Upsilon3Smat[dim3S*dim3S];
	for(int i=0;i<dim2S;i++){
		for(int j=0;j<dim2S;j++){
			int i2=i;
			int j2=j;
			if(j==dim1S) j2=dim2S;
			if(i==dim1S) i2=dim2S;
			Upsilon3Smat[dim2S*i+j]=fullmat[fullmatDim*shift+shift+j2+fullmatDim*i2];
		}
	}

	cout<<"Full matrix"<<endl;
	for(int i=0; i<fullmatDim*fullmatDim; i++){
		printf("\t%d",fullmat[i]);
		if((i+1)%fullmatDim==0 && i>0) printf("\n");
	}
	printf("\n\n");

	cout<<"1S matrix"<<endl;
	for(int i=0; i<dim1S*dim1S; i++){
		printf("\t%d",Upsilon1Smat[i]);
		if((i+1)%dim1S==0 && i>0) printf("\n");
	}
	printf("\n\n");

	cout<<"2S matrix"<<endl;
	for(int i=0; i<dim2S*dim2S; i++){
		printf("\t%d",Upsilon2Smat[i]);
		if((i+1)%dim2S==0 && i>0) printf("\n");
	}
	printf("\n\n");

	cout<<"3S matrix"<<endl;
	for(int i=0; i<dim3S*dim3S; i++){
		printf("\t%d",Upsilon3Smat[i]);
		if((i+1)%dim2S==0 && i>0) printf("\n");
	}
	printf("\n\n");
}