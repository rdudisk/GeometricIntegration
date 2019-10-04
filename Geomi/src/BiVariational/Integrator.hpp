#ifndef DEF_BIVARIATIONAL_INTEGRATOR
#define DEF_BIVARIATIONAL_INTEGRATOR

#include <iostream>
#include <cmath>
#include <vector>
#include <typeinfo> // ??

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
	Abstract::Problem* m_problem;
	//Abstract::Step<T_M,T_Q,T_INTERNALS,T_PROBLEM,T_TQ>& m_step;
	
public:
	Integrator ( )
	{ }

	~Integrator ()
	{ }

	void
	initialize (void)
	{
		// TODO
	}

	bool
	integrate (void)
	{
		bool success = false;

		size_t i,j;
		size_t n_time_steps = m_problem.size(0);
		size_t n_space_steps = m_problem.size(1);

		for (i=1; i<n_time_steps-1; i++) {
			for (j=0; j<n_space_steps-1; j++) {
				// Setting up
				eta_cur[j] = (1.0/s)*T_DIFFEO::inv(m_problem.pos(i,j).inverted()*m_problem->pos(i,j+1));
				lambda_cur[j] = T_DIFFEO::dR_inv_star(s*eta_cur[j],m_problem->dLdv(1,eta_cur[j]));
				mu_cur[j] = 
			}
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
