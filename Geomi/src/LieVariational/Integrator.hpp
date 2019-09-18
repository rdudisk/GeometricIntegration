#ifndef DEF_LIEVARIATIONAL_INTEGRATOR
#define DEF_LIEVARIATIONAL_INTEGRATOR

#include <iostream>
#include <cmath>
#include <vector>
#include <typeinfo> // ??

namespace LieVariational {
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

namespace LieVariational {

template <typename T_M,
		  typename T_Q,
		  typename T_ALGEBRA>
class Integrator : public Abstract::Integrator
{
private:
	Abstract::Problem<T_M,T_Q,T_ALGEBRA>& m_problem;
	Abstract::Step<T_M,T_Q,T_ALGEBRA>& m_step;
	
public:
	Integrator<T_M,T_Q,T_ALGEBRA> (	Abstract::Problem<T_M,T_Q,T_ALGEBRA>& problem,
									Abstract::Step<T_M,T_Q,T_ALGEBRA>& step)
	: m_problem(problem), m_step(step)
	{ }

	~Integrator<T_M,T_Q,T_ALGEBRA> ()
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
