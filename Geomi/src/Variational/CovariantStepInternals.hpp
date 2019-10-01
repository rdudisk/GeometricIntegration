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
	{ return ((1.0/m_h)*T_ALGEBRA::cay_inv(m_q0.inverse()*m_q1)).toNOXVector(); }

	T_Q
	posFromVel (T_M h, T_Q q0, T_ALGEBRA v0) const
	{ return q0*((h*v0).cay()); }

	bool
	computeF (NOXVector<T_ALGEBRA::DOF>& f, const NOXVector<T_ALGEBRA::DOF>& xi)
	{
		T_ALGEBRA xi_prev = (1.0/m_h)*T_ALGEBRA::cay_inv(m_q0.inverse()*m_q1);
		T_ALGEBRA xi_next = T_ALGEBRA(xi);
		T_Q tau_prev = (m_h*xi_prev).cay();

		f =	(-1.0)*T_ALGEBRA::static_Ad_star(tau_prev,T_ALGEBRA((m_h*xi_prev).dCayRInv().transpose()*m_problem.dLdv(xi_prev))).toVector()
				+ (m_h*xi_next).dCayRInv().transpose()*m_problem.dLdv(xi_next);
		return true;
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_ALGEBRA::DOF,T_ALGEBRA::DOF>& J, const NOXVector<T_ALGEBRA::DOF>& xi_vec)
	{
		T_ALGEBRA xi = T_ALGEBRA(xi_vec);

		for (int i=0; i<T_ALGEBRA::DOF; i++) {
			// TODO vérifier les indices
			J.block(0,i,T_ALGEBRA::DOF,1) = (Eigen::Matrix<double,T_ALGEBRA::DOF,T_ALGEBRA::DOF>::Identity()+0.5*m_h*T_ALGEBRA::GeneratorMatrix(i)+m_h*m_h*0.25*(T_ALGEBRA::GeneratorVector(i)*xi.toVector().transpose()+xi.toVector()*T_ALGEBRA::GeneratorVector(i).transpose()))*m_problem.dLdv(xi);
		}
		J += (m_h*xi).dCayRInv().transpose()*m_problem.JvdLdv(xi);

		return true;
	}
};
} // namespace Variational

#endif
