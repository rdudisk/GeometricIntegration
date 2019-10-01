#ifndef DEF_RKMK_EXPLICIT_STEP
#define DEF_RKMK_EXPLICIT_STEP

namespace RKMK {

/**
 * This class inherits RKMK::Abstract::Step and implements an explicit RKMK step.
 * A step is explicit iff the coefficients of the Butcher tableau check \f$ a_{i,j}=0,\,\forall i<j\f$.
 * Since the computation of the solution is straightforward, it is generally faster than an implicit method.
 */

template <typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES,
		  typename T_M = double>
class ExplicitStep : public RKMK::Abstract::Step<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>
{
public:
	ExplicitStep<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M> (Abstract::Problem<T_LIE_ALGEBRA,T_M>& problem)
	:	RKMK::Abstract::Step<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>(problem)
	{ }

	const NOXVector<T_LIE_ALGEBRA::DOF>
	makeStep (void)
	{
		NOXVector<T_LIE_ALGEBRA::DOF>* Y1 = new NOXVector<T_LIE_ALGEBRA::DOF>();
		bool success;

		success = this->m_internals->computeSolution();
		*Y1 = this->m_internals->reconstruct();

		return *Y1;
	}
};
} // namespace RKMK

#endif
