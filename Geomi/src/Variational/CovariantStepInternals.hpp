#ifndef DEF_VARIATIONAL_COVARIANTSTEPINTERNALS
#define DEF_VARIATIONAL_COVARIANTSTEPINTERNALS

#include <vector>

namespace Variational {

template <typename T_M,
		  typename T_Q,
		  typename T_ALGEBRA>
class CovariantStepInternals
	: public Abstract::StepInternals<T_M,T_Q,Abstract::LieProblem<T_M,T_Q,T_ALGEBRA>,T_ALGEBRA>,
	  public ::Abstract::NOXStep<T_ALGEBRA,1>
{
	using Problem = Abstract::LieProblem<T_M,T_Q,T_ALGEBRA>;

protected:
	using Abstract::StepInternals<T_M,T_Q,Problem,T_ALGEBRA>::m_h;
	using Abstract::StepInternals<T_M,T_Q,Problem,T_ALGEBRA>::m_q0;
	using Abstract::StepInternals<T_M,T_Q,Problem,T_ALGEBRA>::m_q1;
	using Abstract::StepInternals<T_M,T_Q,Problem,T_ALGEBRA>::m_problem;

public:
	CovariantStepInternals<T_M,T_Q,T_ALGEBRA> (Problem& problem)
	:	Abstract::StepInternals<T_M,T_Q,Problem,T_ALGEBRA>(problem)
	{ }

	const NOXVector<T_ALGEBRA::DOF>
	getInitialGuess ()
	{
		// TODO
		//NOXVector<T_Q::DOF> ret((1.0+1.0/m_h)*m_q1-(1.0/m_h)*m_q0);
		//return ret;
	}

	T_Q
	posFromVel (T_M h, T_Q q0, T_ALGEBRA v0) const
	{ return q0*((h*v0).cay()); }

	bool
	computeF (NOXVector<T_ALGEBRA::DOF>& f, const NOXVector<T_ALGEBRA::DOF>& q)
	{
		T_ALGEBRA xi_prev = T_ALGEBRA::cay_inv(m_q0.inverse()*m_q1);
		T_ALGEBRA xi_next = T_ALGEBRA(q);
		T_Q tau_prev = xi_prev.cay();
		T_Q tau_next = xi_next.cay();

		T_ALGEBRA::static_bracket(xi_prev,xi_next);
		xi_prev.inverted();
		f =	T_ALGEBRA::static_Ad_star(tau_prev,T_ALGEBRA(xi_prev.dCayRInv().transpose()*m_problem.dLdv(xi_prev*(1.0/m_h)))).toVector()
				- xi_next.dCayRInv().transpose()*m_problem.dLdv(xi_next*(1.0/m_h));
		return true;
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_ALGEBRA::DOF,T_ALGEBRA::DOF>& J, const NOXVector<T_ALGEBRA::DOF>& q)
	{
		//J =		m_problem.JvdLdq(m_q1,(q-m_q1)/m_h )
			//-	m_problem.JvdLdv(m_q1,(q-m_q1)/m_h )/m_h;
		return true;
	}
};
} // namespace Variational

#endif
