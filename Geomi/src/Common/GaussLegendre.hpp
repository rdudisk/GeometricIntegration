#ifndef DEF_INTERPOLATION_GAUSSLEGENDRE
#define DEF_INTERPOLATION_GAUSSLEGENDRE

#include <vector>
#include <cmath>

class GaussLegendre
{
private:
	/* coeffs given for [-1,1] interpolation */
	static const double GL_WEIGHTS[]; 
	static const double GL_DATES[];
	static const int MAX_DEGREE = 5;

public:
	static std::vector<double>
	weights (int degree)
	{
		if (degree>MAX_DEGREE)
			return std::vector<double>();
		std::vector<double> ret;
		// translate from [-1,1] to [0,1] interpolation
		for (int i=0; i<degree; i++) {
			ret.push_back(GL_WEIGHTS[((degree*(degree-1))/2)+i]/2.0);
		}
		return ret;
	}

	static std::vector<double>
	dates (int degree)
	{
		if (degree>MAX_DEGREE)
			return std::vector<double>();
		std::vector<double> ret;
		// translate from [-1,1] to [0,1] interpolation
		for (int i=0; i<degree; i++) {
			ret.push_back((GL_DATES[((degree*(degree-1))/2)+i]/2.0)+0.5);
		}
		return ret;
	}
};

const double GaussLegendre::GL_WEIGHTS[] = {
	2.0,
    1.0,1.0,
    5.0/9.0,8.0/9.0,5.0/9.0,
    (18.0-sqrt(30.0))/36.0,(18.0+sqrt(30.0))/36.0,(18.0+sqrt(30.0))/36.0,(18.0-sqrt(30.0))/36.0,
    (322.0-13.0*sqrt(70.0))/900.0,(322.0+13.0*sqrt(70.0))/900.0,128.0/225.0,(322.0+13.0*sqrt(70.0))/900.0,(322.0-13.0*sqrt(70.0))/900.0 };
const double GaussLegendre::GL_DATES[] = {
	0.0,
	-1.0/sqrt(3),1.0/sqrt(3),
	-sqrt(3.0/5.0),0.0,sqrt(3.0/5.0),
	-sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0)),-sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0)),sqrt((3.0/7.0)-(2.0/7.0)*sqrt(6.0/5.0)),sqrt((3.0/7.0)+(2.0/7.0)*sqrt(6.0/5.0)),
	-(1.0/3.0)*sqrt(5.0+2.0*sqrt(10.0/7.0)),-(1.0/3.0)*sqrt(5.0-2.0*sqrt(10.0/7.0)),0.0,(1.0/3.0)*sqrt(5.0-2.0*sqrt(10.0/7.0)),(1.0/3.0)*sqrt(5.0+2.0*sqrt(10.0/7.0))};

#endif
