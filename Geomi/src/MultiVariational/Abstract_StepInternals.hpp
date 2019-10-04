#ifndef DEF_MULTIVARIATIONAL_ABSTRACT_STEPINTERNALS
#define DEF_MULTIVARIATIONAL_ABSTRACT_STEPINTERNALS

#include <vector>

namespace MultiVariational {
namespace Abstract {

template <typename T_Q,
		  typename T_PROBLEM,
		  int T_N_INDEP,
		  typename T_TQ = T_Q>
class StepInternals
{
protected:
	std::vector<double> m_h;
	T_Q m_q0;
	std::vector<T_Q> m_q1;

	T_PROBLEM& m_problem;

public:
	StepInternals<T_Q,T_PROBLEM,T_N_INDEP,T_TQ> (T_PROBLEM& problem)
	: m_problem(problem)
	{ 
		m_h = std::vector<double>(T_N_INDEP);
		m_q1 = std::vector<T_Q>(T_N_INDEP);
	}

	void
	setData (std::vector<double> h, T_Q q0, std::vector<T_Q> q1)
	{
		// TODO: check sizes
		m_h = h;
		m_q0 = q0;
		m_q1 = q1;
	}

	std::vector<double>
	h () const
	{ return m_h; }

	T_Q
	q0 () const
	{ return m_q0; }

	std::vector<T_Q>
	q1 () const
	{ return m_q1; }

	virtual T_Q
	posFromVel (double h, T_Q q0, T_TQ v0) const = 0;
};
} // namespace Abstract
} // namespace MultiVariational

#endif
