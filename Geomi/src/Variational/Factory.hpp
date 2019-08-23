#ifndef DEF_VARIATIONAL_FACTORY
#define DEF_VARIATIONAL_FACTORY
	
namespace Variational {

/**
 * %Factory for Variational integrators.
 * The following integrators are implemented:
 *
 *	- Explicit Euler
 */

template <typename T_M,
		  typename T_Q,
		  typename T_TQ>
class Factory
{
public:
	static Abstract::Integrator&
	createIntegrator (Abstract::Problem<T_M,T_Q>& problem, std::string s)
	{
		Abstract::Integrator* integrator = NULL;

		if (s=="Explicit Euler") {
			EulerStep<T_M,T_Q,T_TQ>* step = new EulerStep<T_M,T_Q,T_TQ>(problem);
			integrator = new Integrator<T_M,T_Q,T_TQ>(problem, *step);
		}
		else if (s=="Implicit Euler") {
			ImplicitEulerStep<T_M,T_Q,T_TQ>* step = new ImplicitEulerStep<T_M,T_Q,T_TQ>(problem);
			integrator = new Integrator<T_M,T_Q,T_TQ>(problem, *step);
		}
		else if (s=="Midpoint") {
			MidpointStep<T_M,T_Q,T_TQ>* step = new MidpointStep<T_M,T_Q,T_TQ>(problem);
			integrator = new Integrator<T_M,T_Q,T_TQ>(problem, *step);
		}
		// notation PsNrGau
		// vient de PsNrQu:
		// polynomials of degree s with s+1 control points, quadrature formula of order u with r quadrature points (with gauss u=2r)
		else if (s=="Galerkin P2N2Gau") {
			GalerkinStep<T_M,T_Q,T_TQ,2>* step = new GalerkinStep<T_M,T_Q,T_TQ,2>(problem);
			integrator = new Integrator<T_M,T_Q,T_TQ>(problem, *step);
		}

		return *integrator;
	}
};

} // namespace Variational

#endif

