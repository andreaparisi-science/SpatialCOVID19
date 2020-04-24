#include "SplineInterpolator.h"


SplineInterpolator::SplineInterpolator( const std::vector<double> &xx, const std::vector<double> &yy, const int method, const std::vector<double> &yyfd )  {
	int nn = xx.size();
	double sig, pp, qqn, uun;
	std::vector<double>  uu(nn-1, 0);

	yysd.resize(nn, 0);
	if (method == SPLINE_METHOD_NATURAL)  {
		yysd[0] = uu[0] = 0.0;
	} else {
		yysd[0] = -0.5;
		uu[0] = (3.0/(xx[1]-xx[0]))*((yy[1]-yy[0])/(xx[1]-xx[0])-yyfd[0]);
	}
	for (int jj = 1; jj < nn; jj++)  {
		sig = (xx[jj]-xx[jj-1])/(xx[jj+1]-xx[jj-1]);
		pp = sig*yysd[jj-1]+2.0;
		yysd[jj] = (sig-1.0)/pp;
		uu[jj] = (yy[jj+1]-yy[jj])/(xx[jj+1]-xx[jj]) - (yy[jj]-yy[jj-1])/(xx[jj]-xx[jj-1]);
		uu[jj] = (6.0*uu[jj]/(xx[jj+1]-xx[jj-1])-sig*uu[jj-1])/pp;
	}
	if (method == SPLINE_METHOD_NATURAL)  {
		qqn = uun = 0.0;
	} else {
		qqn = +0.5;
		uun = (3.0/(xx[nn-1]-xx[nn-2]))*(yyfd[1]-(yy[nn-1]-yy[nn-2])/(xx[nn-1]-xx[nn-2]));
	}
	yysd[nn-1] = (uun-qqn*uu[nn-2])/(qqn*yysd[nn-2]+1.0);
	for (int kk=nn-2; kk>=0; kk--)  {
		yysd[kk]=yysd[kk]*yysd[kk+1]+uu[kk];
	}
	xval = xx;
	yval = yy;
}



SplineInterpolator::~SplineInterpolator()  {
}



double  SplineInterpolator::interpolate( double xx )  {
	int nn = xval.size();
	int kk, kklo = 0;
	int kkhi = nn-1;
	double hh, aa, bb;
	while ( kkhi-kklo > 1 )  {
		kk = (kkhi+kklo) >> 1;
		if (xval[kk] > xx)  kkhi = kk;
		else kklo = kk;
	}
	hh = xval[kkhi] - xval[kklo];
	if (hh == 0.0)  {
		throw  "Bad set of x values.";
	}
	aa = (xval[kkhi]-xx)/hh;
	bb = (xx-xval[kklo])/hh;
	return  aa*yval[kklo]+bb*yval[kkhi]+ ((aa*aa*aa-aa)*yysd[kklo]+(bb*bb*bb-bb)*yysd[kkhi])*(hh*hh)/6.0;
}



std::vector<double>  SplineInterpolator::interpolate( const std::vector<double> &xx )  {
	std::vector<double>  res( xx.size(), 0 );
	for (int jj = 0; jj < xx.size(); jj++)  {
		res[jj] = interpolate( xx[jj] );
	}
	return res;
}


