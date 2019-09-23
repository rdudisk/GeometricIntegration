#ifndef DEF_VARIATIONAL_ABSTRACT_STEP
#define DEF_VARIATIONAL_ABSTRACT_STEP

namespace Variational {
namespace Abstract {

template <typename T_M,
		  typename T_Q,
		  typename T_INTERNALS,
		  typename T_PROBLEM,
		  typename T_TQ = T_Q>
class Step
{ 
protected:
	T_INTERNALS* m_internals;
	T_PROBLEM& m_problem;

public:
	Step<T_M,T_Q,T_INTERNALS,T_PROBLEM> (T_PROBLEM& problem)
	:	m_problem(problem)
	{ m_internals = new T_INTERNALS(problem); }

	~Step<T_M,T_Q,T_INTERNALS,T_PROBLEM> ()
	{ }

	T_Q
	posFromVel (T_Q q0, T_TQ v0)
	{ return m_internals->posFromVel(q0,v0); }

	virtual void
	setData (T_M h_var, T_Q q0_var, T_Q q1_var) = 0;

	virtual const T_Q
	makeStep () = 0;

	virtual void
	initialize () = 0;
};
} // namespace Abstract
} // namespace Variational

#endif
