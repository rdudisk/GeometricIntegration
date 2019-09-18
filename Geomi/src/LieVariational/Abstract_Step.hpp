#ifndef DEF_LIEVARIATIONAL_ABSTRACT_STEP
#define DEF_LIEVARIATIONAL_ABSTRACT_STEP

namespace LieVariational {
namespace Abstract {

template <typename T_M,
		  typename T_Q,
		  typename T_ALGEBRA>
class Step
{ 
protected:
	Variational::Abstract::StepInternals<T_M,T_Q,T_ALGEBRA>* m_internals;
	Variational::Abstract::Problem<T_M,T_Q,T_ALGEBRA>& m_problem;

public:
	Step<T_M,T_Q,T_ALGEBRA> (Variational::Abstract::Problem<T_M,T_ALGEBRA>& problem)
	:	m_problem(problem)
	{ m_internals = new Variational::Abstract::StepInternals<T_M,T_Q,T_ALGEBRA>(problem); }

	~Step<T_M,T_Q,T_ALGEBRA> ()
	{ }

	T_Q
	posFromVel (T_Q q0, T_ALGEBRA v0)
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
