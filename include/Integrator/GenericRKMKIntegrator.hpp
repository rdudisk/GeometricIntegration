#ifndef DEF_GENERIC_RKMK_INTEGRATOR
#define DEF_GENERIC_RKMK_INTEGRATOR

#include "include/Problem/ProblemInterface.hpp"
#include "include/Step/GenericRKMKStep.hpp"
#include "include/Common/Utils.hpp"

namespace Abstract {
class RKMKIntegrator
{
public:
	virtual void
	setCoeffs (std::vector<double> va, std::vector<double> vb) = 0;

	virtual bool
	integrate () = 0;
};
} // namespace Abstract

template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES>
class RKMKIntegrator : public Abstract::RKMKIntegrator
{
private:
	Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& m_problem;
	Abstract::RKMKStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>& m_step;
	
public:
	RKMKIntegrator<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& problem,
															   Abstract::RKMKStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>& step)
	: m_problem(problem), m_step(step)
	{ }

	~RKMKIntegrator<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> ()
	{ }

	void
	setCoeffs (std::vector<double> va, std::vector<double> vb)
	{ m_step.setCoeffs(va,vb); }

	bool
	integrate (void)
	{
		int i;
		bool success = false;
		double h;
		int n_steps = m_problem.size();
		T_Q Y0, Y1;

		for (i=0; i<n_steps-1; i++) {
			Y0 = m_problem.pos(i);
			h = m_problem.base(i+1)-m_problem.base(i);
			m_step.setData(h,Y0);

			Y1 = m_step.makeStep();
			m_problem.pos(i+1,Y1);
		}

		return success;
	}
}; 
	
template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA>
class RKMKFactory
{
public:
	static Abstract::RKMKIntegrator&
	createIntegrator (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& problem, std::string s)
	{
		Abstract::RKMKIntegrator* integrator = NULL;

		if (s=="Explicit Euler") {
			double a[] = { 0.0 };
			double b[] = { 1.0 };
			std::vector<double> va(a,a+sizeof(a)/sizeof(double));
			std::vector<double> vb(b,b+sizeof(b)/sizeof(double));
			ExplicitRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,1>* step = new ExplicitRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,1>(problem);
			step->setCoeffs(va,vb);
			integrator = new RKMKIntegrator<T_M,T_Q,T_LIE_ALGEBRA,1>(problem, *step);
			integrator->setCoeffs(va,vb);
		}
		else if (s=="Implicit Euler") {
			double a[] = { 1.0 };
			double b[] = { 1.0 };
			std::vector<double> va(a,a+sizeof(a)/sizeof(double));
			std::vector<double> vb(b,b+sizeof(b)/sizeof(double));
			DIRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,1>* step = new DIRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,1>(problem);
			step->setCoeffs(va,vb);
			integrator = new RKMKIntegrator<T_M,T_Q,T_LIE_ALGEBRA,1>(problem, *step);
			integrator->setCoeffs(va,vb);
		}
		else if (s=="Midpoint") {
			double a[] = { 0.5 };
			double b[] = { 1.0 };
			std::vector<double> va(a,a+sizeof(a)/sizeof(double));
			std::vector<double> vb(b,b+sizeof(b)/sizeof(double));
			DIRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,1>* step = new DIRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,1>(problem);
			step->setCoeffs(va,vb);
			integrator = new RKMKIntegrator<T_M,T_Q,T_LIE_ALGEBRA,1>(problem, *step);
			integrator->setCoeffs(va,vb);
		}
		else if (s=="RK 4") {
			double a[] = { 0.0,0.0,0.0,0.0, 0.5,0.0,0.0,0.0, 0.0,0.5,0.0,0.0, 0.0,0.0,1.0,0.0 };
			double b[] = { 1.0/6.0,2.0/6.0,2.0/6.0,1.0/6.0 };
			std::vector<double> va(a,a+sizeof(a)/sizeof(double));
			std::vector<double> vb(b,b+sizeof(b)/sizeof(double));
			ExplicitRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,4>* step = new ExplicitRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,4>(problem);
			step->setCoeffs(va,vb);
			integrator = new RKMKIntegrator<T_M,T_Q,T_LIE_ALGEBRA,4>(problem, *step);
			integrator->setCoeffs(va,vb);
		}

		return *integrator;
	}
};

#endif
