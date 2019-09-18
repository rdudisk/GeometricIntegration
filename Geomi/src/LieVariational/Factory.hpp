#ifndef DEF_LIEVARIATIONAL_FACTORY
#define DEF_LIEVARIATIONAL_FACTORY
	
namespace LieVariational {

template <typename T_M,
		  typename T_Q,
		  typename T_ALGEBRA>
class Factory
{
public:
	static Abstract::Integrator&
	createIntegrator (Abstract::Problem<T_M,T_Q,T_ALGEBRA>& problem, std::string s)
	{
		Abstract::Integrator* integrator = NULL;

		if (s=="Euler") {
			EulerStep<T_M,T_Q,T_TQ>* step
				= new EulerStep<T_M,T_Q,T_TQ>(problem);
			integrator = new Integrator<T_M,T_Q,T_TQ>(problem, *step);
		}
		return *integrator;
	}
};

} // namespace Variational

#endif

