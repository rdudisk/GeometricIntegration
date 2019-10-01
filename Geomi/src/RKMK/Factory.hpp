#ifndef DEF_RKMK_FACTORY
#define DEF_RKMK_FACTORY

//#include "RKMK_Abstract_Integrator.hpp"
//#include "RKMK_ExplicitStep.hpp"
//#include "RKMK_DiagonalStep.hpp"
	
namespace RKMK {

/**
 * %Factory for RKMK integrators.
 * This is meant to be a user friendly interface to instantiate RKMK::Integrator objects.
 * The following integrators are implemented:
 *
 *	- Explicit Euler
 *	- Implicit Euler
 *	- Midpoint
 *	- Classical RK 4
 *
 * If you wish to use an integrator that is not listed, please refer to the documentation for RKMK::Integrator.
 */

template <typename T_LIE_ALGEBRA, typename T_M = double>
class Factory
{
public:
	static Abstract::Integrator&
	createIntegrator (Abstract::Problem<T_LIE_ALGEBRA,T_M>& problem, std::string s)
	{
		Abstract::Integrator* integrator = NULL;

		if (s=="Explicit Euler") {
			double a[] = { 0.0 };
			double b[] = { 1.0 };
			std::vector<double> va(a,a+sizeof(a)/sizeof(double));
			std::vector<double> vb(b,b+sizeof(b)/sizeof(double));
			ExplicitStep<T_LIE_ALGEBRA,1,T_M>* step = new ExplicitStep<T_LIE_ALGEBRA,1,T_M>(problem);
			integrator = new Integrator<T_LIE_ALGEBRA,1,T_M>(problem, *step);
			integrator->setCoeffs(va,vb);
		}
		else if (s=="Implicit Euler") {
			double a[] = { 1.0 };
			double b[] = { 1.0 };
			std::vector<double> va(a,a+sizeof(a)/sizeof(double));
			std::vector<double> vb(b,b+sizeof(b)/sizeof(double));
			DiagonalStep<T_LIE_ALGEBRA,1,T_M>* step = new DiagonalStep<T_LIE_ALGEBRA,1,T_M>(problem);
			integrator = new Integrator<T_LIE_ALGEBRA,1,T_M>(problem, *step);
			integrator->setCoeffs(va,vb);
		}
		else if (s=="Midpoint") {
			double a[] = { 0.5 };
			double b[] = { 1.0 };
			std::vector<double> va(a,a+sizeof(a)/sizeof(double));
			std::vector<double> vb(b,b+sizeof(b)/sizeof(double));
			DiagonalStep<T_LIE_ALGEBRA,1,T_M>* step = new DiagonalStep<T_LIE_ALGEBRA,1,T_M>(problem);
			integrator = new Integrator<T_LIE_ALGEBRA,1,T_M>(problem, *step);
			integrator->setCoeffs(va,vb);
		}
		else if (s=="RK 4") {
			double a[] = { 0.0,0.0,0.0,0.0, 0.5,0.0,0.0,0.0, 0.0,0.5,0.0,0.0, 0.0,0.0,1.0,0.0 };
			double b[] = { 1.0/6.0,2.0/6.0,2.0/6.0,1.0/6.0 };
			std::vector<double> va(a,a+sizeof(a)/sizeof(double));
			std::vector<double> vb(b,b+sizeof(b)/sizeof(double));
			ExplicitStep<T_LIE_ALGEBRA,4,T_M>* step = new ExplicitStep<T_LIE_ALGEBRA,4,T_M>(problem);
			integrator = new Integrator<T_LIE_ALGEBRA,4,T_M>(problem, *step);
			integrator->setCoeffs(va,vb);
		}

		return *integrator;
	}
};

} // namespace RKMK

#endif
