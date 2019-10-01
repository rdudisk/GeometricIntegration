#ifndef DEF_RKMK_ABSTRACT_INTEGRATOR
#define DEF_RKMK_ABSTRACT_INTEGRATOR

namespace RKMK {
namespace Abstract {

/**
 * Abstract class for RKMK integrators.
 */

class Integrator
{
public:
	virtual void
	setCoeffs (std::vector<double> va, std::vector<double> vb) = 0;

	virtual bool
	integrate () = 0;
};
} // namespace Abstract
} // namespace RKMK

namespace RKMK {

/**
 * Template class for RKMK integrators.
 * This class should not be instatiated explicitely by the user, but used through RKMK::Factory instead.
 * If RKMK::Factory does not provide the appropriate integrator, here is a code snippet for the classical RK 4 method to help you create your own.
 *	\code{.cpp}
	// Instance of class MyProblem implementing Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>
	MyProblem problem;

	// Butcher tableau coefficients
	double a[] = { 0.0,0.0,0.0,0.0, 0.5,0.0,0.0,0.0, 0.0,0.5,0.0,0.0, 0.0,0.0,1.0,0.0 };
	double b[] = { 1.0/6.0,2.0/6.0,2.0/6.0,1.0/6.0 };
	std::vector<double> va(a,a+sizeof(a)/sizeof(double));
	std::vector<double> vb(b,b+sizeof(b)/sizeof(double));

	// Define the appropriate step
	RKMK::ExplicitStep<T_M,T_Q,T_LIE_ALGEBRA,4> step(problem);

	// Create the integrator form the problem and the step
	RKMK::Integrator<T_M,T_Q,T_LIE_ALGEBRA,4> integrator(problem, step);
	integrator.setCoeffs(va,vb);
	\endcode
 */

template <typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES,
		  typename T_M = double>
class Integrator : public RKMK::Abstract::Integrator
{
private:
	RKMK::Abstract::Problem<T_LIE_ALGEBRA,T_M>& m_problem;
	RKMK::Abstract::Step<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>& m_step;
	
public:
	Integrator<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M> (RKMK::Abstract::Problem<T_LIE_ALGEBRA,T_M>& problem,
													   RKMK::Abstract::Step<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>& step)
	: m_problem(problem), m_step(step)
	{ }

	~Integrator<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M> ()
	{ }

	void
	setCoeffs (std::vector<double> va, std::vector<double> vb)
	{ m_step.setCoeffs(va,vb); }

	bool
	integrate (void)
	{
		int i;
		bool success = false;
		T_M h;
		int n_steps = m_problem.size();
		NOXVector<T_LIE_ALGEBRA::DOF> Y0, Y1;

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
} // namespace RKMK

#endif
