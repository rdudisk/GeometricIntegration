#ifndef DEF_BIVARIATIONAL_INTEGRATOR
#define DEF_BIVARIATIONAL_INTEGRATOR

#include <iostream>
#include <cmath>
#include <vector>
#include <typeinfo>

/*
namespace BiVariational {
namespace Abstract {

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
*/

namespace BiVariational {

template <typename T_SCALAR,
		  typename T_Q,
		  typename T_VEL>
		  //typename T_INTERNALS,
		  //typename T_PROBLEM>
class Integrator
{
private:
	Abstract::LieProblem<T_SCALAR,T_Q,T_VEL>* m_problem;
	Epetra_Comm* m_comm;
	Epetra_Map* m_map;
	Step<T_SCALAR,T_Q,T_VEL>* m_step;
	
public:
	Integrator (Abstract::LieProblem<T_SCALAR,T_Q,T_VEL>& problem, Epetra_Comm& c)
	:	m_problem(&problem),
		m_comm(&c)
	{
		m_step = new Step<T_SCALAR,T_Q,T_VEL>(*m_problem,*m_comm);
	}

	~Integrator ()
	{ }

	bool
	initialize (void)
	{
		bool success = m_step->initialize();
		return success;
	}

	bool
	integrate (void)
	{
		bool success = true;

		//m_step->test();

		
		size_t i,j;
		size_t n_time_steps = m_problem->size(0);
		size_t n_space_steps = m_problem->size(1);

		/* Epetra initialisations */
		const int numGlobEntries = T_Q::DOF * n_space_steps;
		const int indexBase = 0;

		for (i=1; i<n_time_steps-1; i++) {
			m_step->makeStep();
		}

		return success;
	}
}; 
} // namespace Variational

#endif
