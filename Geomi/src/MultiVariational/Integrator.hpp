#ifndef DEF_MULTIVARIATIONAL_INTEGRATOR
#define DEF_MULTIVARIATIONAL_INTEGRATOR

#include <iostream>
#include <cmath>
#include <vector>
#include <typeinfo> // ??

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


namespace MultiVariational {
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
} // namespace MultiVariational

namespace MultiVariational {

template <typename T_Q,
		  typename T_INTERNALS,
		  typename T_PROBLEM,
		  int T_N_INDEP,
		  typename T_TQ = T_Q>
class Integrator : public Abstract::Integrator
{
private:
	T_PROBLEM& m_problem;
	Abstract::Step<T_Q,T_INTERNALS,T_PROBLEM,T_N_INDEP,T_TQ>& m_step;
	std::vector<int> m_direction;
	
public:
	Integrator<T_Q,T_INTERNALS,T_PROBLEM,T_N_INDEP,T_TQ> (T_PROBLEM& problem,
		  Abstract::Step<T_Q,T_INTERNALS,T_PROBLEM,T_N_INDEP,T_TQ>& step)
	: m_problem(problem), m_step(step)
	{ }

	~Integrator<T_Q,T_INTERNALS,T_PROBLEM,T_TQ> ()
	{ }

	void
	setDirection (std::vector<int> direction)
	{
		// TODO: checks
		m_direction = direction;
	}

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
		// Pour l'instant seulement space, then time integration
		int i,j;
		int index[2];
		bool success = false;
		std::vector<double> coord;
		std::vector<int> n_steps = m_problem.size();
		std::vector<MultiVariational::Syst<T_Q,T_N_INDEP>> node_list;

		for (i=1; i<n_steps[0]-1; i++) {
			for (j=1; j<n_steps[1]-1; j++) {
				node_list.clear();
				node_list.push_back(m_problem.node({i,j}));
				node_list.push_back(m_problem.node({i,j}));
				node_list.push_back(m_problem.node({i,j}));
				node_list.push_back(m_problem.node({i,j}));
				node_list.push_back(m_problem.node({i,j}));
				m_step.setData(node_list);

				// ! le schema doit etre un space then time aussi
				q1 = m_step.makeStep();
				m_problem.pos(i+2,q1);
			}
		}

		return success;
	}
}; 
} // namespace MultiVariational

#endif
