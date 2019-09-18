#ifndef DEF_LIEVARIATIONAL_ABSTRACT_STEPINTERNALS
#define DEF_LIEVARIATIONAL_ABSTRACT_STEPINTERNALS

#include <vector>

namespace LieVariational {
namespace Abstract {

template <typename T_M,
		  typename T_Q,
		  typename T_ALGEBRA>
class StepInternals
{
protected:
	T_M m_h;
	T_Q m_q0;
	T_Q m_q1;

	Abstract::Problem<T_M,T_Q,T_ALGEBRA>& m_problem;

public:
	StepInternals<T_M,T_Q,T_ALGEBRA> (Abstract::Problem<T_M,T_Q,T_ALGEBRA>& problem)
	: m_problem(problem)
	{ }

	void
	setData (T_M h, T_Q q0, T_Q q1)
	{
		m_h = h;
		m_q0 = q0;
		m_q1 = q1;
	}
};
} // namespace Abstract
} // namespace Variational

#endif
