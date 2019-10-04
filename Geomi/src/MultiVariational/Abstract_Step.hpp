#ifndef DEF_MULTIVARIATIONAL_ABSTRACT_STEP
#define DEF_MULTIVARIATIONAL_ABSTRACT_STEP

namespace MultiVariational {
namespace Abstract {

template <typename T_Q,
		  typename T_INTERNALS,
		  typename T_PROBLEM,
		  typename T_TQ = T_Q>
class Step
{ 
protected:
	T_INTERNALS* m_internals;
	T_PROBLEM& m_problem;

public:
	Step<T_Q,T_INTERNALS,T_PROBLEM,T_TQ> (T_PROBLEM& problem)
	:	m_problem(problem)
	{ m_internals = new T_INTERNALS(problem); }

	~Step<T_Q,T_INTERNALS,T_PROBLEM,T_TQ> ()
	{ }

	T_Q
	posFromVel (double h, T_Q q0, T_TQ v0)
	{ return m_internals->posFromVel(h,q0,v0); }

	virtual void
	setData (std::vector<double> h_var, T_Q q0_var, std::vector<T_Q> q1_var) = 0;

	virtual const T_Q
	makeStep () = 0;

	virtual void
	initialize () = 0;
};
} // namespace Abstract
} // namespace MultiVariational

#endif
