#ifndef DEF_VARIATIONAL_ABSTRACT_STEP
#define DEF_VARIATIONAL_ABSTRACT_STEP

namespace Variational {
namespace Abstract {

template <typename T_M,
		  typename T_Q,
		  typename T_TQ>
class Step
{ 
protected:
	Variational::Abstract::StepInternals<T_M,T_Q,T_TQ>* m_internals;
	Variational::Abstract::Problem<T_M,T_Q>& m_problem;

public:
	Step<T_M,T_Q,T_TQ> (Variational::Abstract::Problem<T_M,T_Q>& problem)
	:	m_problem(problem)
	{
		m_internals = new Variational::Abstract::StepInternals<T_M,T_Q,T_TQ>(problem);
		//setType();
	}

	~Step<T_M,T_Q,T_TQ> ()
	{ }

	T_Q
	posFromVel (T_Q q0, T_TQ v0)
	{
		return m_internals->posFromVel(q0,v0);
	}

	virtual void
	setData (T_M h_var, T_Q q0_var, T_Q q1_var) = 0;
	//{ m_internals->setData(h_var,q0_var,q1_var); }

	virtual const T_Q
	makeStep () = 0;

	virtual void
	initialize () = 0;
};
} // namespace Abstract
} // namespace Variational

#endif
