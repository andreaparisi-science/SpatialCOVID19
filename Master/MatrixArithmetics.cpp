void LUdecomposition(std::vector< std::vector<double> > &mm, int nn, std::vector<int> &indx, double *d)  {
	int    imax;
	double big,dum,sum,temp;
	double std::vector<double> vect(nn, 0);

	*d=1.0;
	for (int ii=0; ii < nn; ii++) {
		big=0.0;
		for (int jj=0; jj < nn; jj++)  {
			if ((temp=fabs(mm[ii][jj])) > big) big=temp;
		}
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vect[ii]=1.0/big;
	}
	for (int jj=0; jj < nn; jj++) {
		for (int ii=0; ii<jj; ii++) {
			sum=mm[ii][jj];
			for (int kk=0; kk<ii; kk++) sum -= a[ii][kk]*a[kk][jj];
			a[ii][jj]=sum;
		}
		big=0.0;
		for (int ii=jj; ii < nn; ii++) {
			sum=mm[ii][jj];
			for (int kk=0; kk<jj; kk++)  sum -= mm[ii][kk]*mm[kk][jj];
			mm[ii][jj] = sum;
			if ( (dum=vect[ii]*fabs(sum)) >= big) {
				big=dum;
				imax=ii;
			}
		}
		if (jj != imax) {
			for (int kk=0; kk < nn; kk++) {
				dum=mm[imax][kk];
				mm[imax][kk]=mm[jj][kk];
				mm[jj][kk]=dum;
			}
			*d = -(*d);
			vect[imax]=vect[jj];
		}
		indx[jj]=imax;
		if (mm[jj][jj] == 0.0) mm[jj][jj]=TINY;
		if (jj != nn) {
			dum=1.0/(mm[jj][jj]);
			for (ii=jj+1; ii<nn; ii++)  {
				mm[ii][jj] *= dum;
			}
		}
	}
}



void LUsubstitute(float **a, int n, int *indx, float b[])  {
	int i,ii=0,ip,j;
	float sum;
	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


// Inverse
#define N ...
float **a,**y,d,*col;
int i,j,*indx;
...
ludcmp(a,N,indx,&d);
for(j=1;j<=N;j++) {
	for(i=1;i<=N;i++) col[i]=0.0;
	col[j]=1.0;
	lubksb(a,N,indx,col);
	for(i=1;i<=N;i++) y[i][j]=col[i];
}

// Determinant
#define N ...
float **a,d;
int j,*indx;
...
ludcmp(a,N,indx,&d);
for(j=1;j<=N;j++) d *= a[j][j];

