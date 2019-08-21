#ifndef DEF_RKMK_EXPLICIT_STEP
#define DEF_RKMK_EXPLICIT_STEP

//#include "RKMK_Abstract_Step.hpp"

namespace RKMK {

/**
 * This class inherits RKMK::Abstract::Step and implements an explicit RKMK step.
 * A step is explicit iff the coefficients of the Butcher tableau check \f$ a_{i,j}=0,\,\forall i<j\f$.
 * Since the computation of the solution is straightforward, it is generally faster than an implicit method.
 */

template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES>
class ExplicitStep : public RKMK::Abstract::Step<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>
{
public:
	ExplicitStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface)
	:	RKMK::Abstract::Step<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>(interface)
	{ }

	const T_Q
	makeStep (void)
	{
		T_Q* Y1 = new T_Q();
		bool success;

		success = this->m_internals->computeSolution();
		*Y1 = this->m_internals->reconstruct();

		return *Y1;
	}
};
} // namespace RKMK

#endif
