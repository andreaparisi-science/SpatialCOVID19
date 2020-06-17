#ifndef  H_LSE_SPLINEINTERPOLATOR
#define  H_LSE_SPLINEINTERPOLATOR

#include <vector>

int constexpr  SPLINE_METHOD_NATURAL = 1;
int constexpr  SPLINE_METHOD_FIXDERIV = 2;

class  SplineInterpolator  {
public:
	SplineInterpolator( const std::vector<double> &xx, const std::vector<double> &yy, const int method, const std::vector<double> &yyfd = {0, 0} );
	~SplineInterpolator();
	std::vector<double>  interpolate( const std::vector<double> &xx );
	double  interpolate( double xx );
private:
	std::vector<double>  xval;
	std::vector<double>  yval;
	std::vector<double>  yysd;
};

#endif



