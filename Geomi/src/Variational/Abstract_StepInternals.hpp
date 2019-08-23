#ifndef DEF_VARIATIONAL_ABSTRACT_STEPINTERNALS
#define DEF_VARIATIONAL_ABSTRACT_STEPINTERNALS

#include <vector>

namespace Variational {
namespace Abstract {

template <typename T_M,
		  typename T_Q,
		  typename T_TQ>
class StepInternals
{
protected:
	T_M m_h;
	T_Q m_q0;
	T_Q m_q1;

	Abstract::Problem<T_M,T_Q>& m_problem;

public:
	StepInternals<T_M,T_Q,T_TQ> (Abstract::Problem<T_M,T_Q>& problem)
	:	m_problem(problem)
	{ }

	void
	setData (T_M h, T_Q q0, T_Q q1)
	{
		m_h = h;
		m_q0 = q0;
		m_q1 = q1;
	}

	//virtual T_Q
	//posFromVel (T_Q, T_TQ v0) = 0;
};
} // namespace Abstract
} // namespace Variational

#endif
