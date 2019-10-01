#ifndef DEF_VARIATIONAL_INTEGRATOR
#define DEF_VARIATIONAL_INTEGRATOR

#include <iostream>
#include <cmath>
#include <vector>
#include <typeinfo> // ??

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


namespace Variational {
namespace Abstract {

/**
 * Abstract class for Variational integrators.
 */

class Integrator
{
public:
	virtual bool
	integrate () = 0;

	virtual void
	initialize () = 0;
};
} // namespace Abstract
} // namespace Variational

namespace Variational {

template <typename T_M,
		  typename T_Q,
		  typename T_INTERNALS,
		  typename T_PROBLEM,
		  typename T_TQ = T_Q>
class Integrator : public Abstract::Integrator
{
private:
	T_PROBLEM& m_problem;
	Abstract::Step<T_M,T_Q,T_INTERNALS,T_PROBLEM,T_TQ>& m_step;
	
public:
	Integrator<T_M,T_Q,T_INTERNALS,T_PROBLEM,T_TQ> (T_PROBLEM& problem,
		  Abstract::Step<T_M,T_Q,T_INTERNALS,T_PROBLEM,T_TQ>& step)
	: m_problem(problem), m_step(step)
	{ }

	~Integrator<T_M,T_Q,T_INTERNALS,T_PROBLEM,T_TQ> ()
	{ }

	// is supposed to do nothing for 1 step methods
	void
	initialize (void)
	{
		/* Achtung ! On suppose que pos(0) et pos(1) sont initialisés */
		int i;
		bool success = false;

		T_Q q0 = m_problem.pos(0);
		T_Q q1 = m_problem.pos(1);
		double h = m_problem.base(1)-m_problem.base(0);

		m_step.setData(h,q0,q1);

		m_step.initialize();
	}

	bool
	integrate (void)
	{
		/* Achtung ! On suppose que pos(0) et pos(1) sont initialisés */
		int i;
		bool success = false;
		double h;
		int n_steps = m_problem.size();
		T_Q q0, q1;

		for (i=0; i<n_steps-2; i++) {
			q0 = m_problem.pos(i);
			q1 = m_problem.pos(i+1);
			h = m_problem.base(i+1)-m_problem.base(i);
			m_step.setData(h,q0,q1);

			q1 = m_step.makeStep();
			m_problem.pos(i+2,q1);
		}

		return success;
	}
}; 
} // namespace Variational

#endif
