#ifndef DEF_VARIATIONAL_ABSTRACT_STEPINTERNALS
#define DEF_VARIATIONAL_ABSTRACT_STEPINTERNALS

#include <vector>

namespace Variational {
namespace Abstract {

template <typename T_M,
		  typename T_Q,
		  typename T_PROBLEM>
class StepInternals
{
protected:
	T_M m_h;
	T_Q m_q0;
	T_Q m_q1;

	T_PROBLEM& m_problem;

public:
	StepInternals<T_M,T_Q,T_PROBLEM> (T_PROBLEM& problem)
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
