#ifndef DEF_VARIATIONAL_COVARIANTSTEPINTERNALS
#define DEF_VARIATIONAL_COVARIANTSTEPINTERNALS

#include <vector>

namespace Variational {

template <typename T_M,
		  typename T_Q,
		  typename T_ALGEBRA>
class CovariantStepInternals : public Abstract::StepInternals<T_M,T_Q,Abstract::LieProblem<T_M,T_Q,T_ALGEBRA>>, public ::Abstract::NOXStep<T_ALGEBRA,1>
{
	using Problem = Abstract::LieProblem<T_M,T_Q,T_ALGEBRA>;

public:
	CovariantStepInternals<T_M,T_Q,T_ALGEBRA> (Problem& problem)
	:	Abstract::StepInternals<T_M,T_Q,Problem>(problem)
	{ }

	const NOXVector<T_ALGEBRA::DOF>
	getInitialGuess ()
	{
		NOXVector<T_Q::DOF> ret((1.0+1.0/this->m_h)*this->m_q1-(1.0/this->m_h)*this->m_q0);
		return ret;
	}

	bool
	computeF (NOXVector<T_ALGEBRA::DOF>& f, const NOXVector<T_ALGEBRA::DOF>& q)
	{
		T_ALGEBRA xi_prev = T_ALGEBRA::cay_inv(this->m_q0.inverse()*this->m_q1);
		T_ALGEBRA xi_next = T_ALGEBRA(q);
		T_Q tau_prev = T_ALGEBRA::cay(xi_prev);
		T_Q tau_next = T_ALGEBRA::cay(xi_next);

		f =	T_ALGEBRA::Ad(tau_prev).transpose()*xi_prev.dCayRInv().transpose()*this->m_problem.dLdv(xi_prev/this->m_h)
				- xi_next.dCayRInv().transpose()*this->m_problem.dLdv(xi_next/this->m_h);
		return true;
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_ALGEBRA::DOF,T_ALGEBRA::DOF>& J, const NOXVector<T_ALGEBRA::DOF>& q)
	{
		//J =		this->m_problem.JvdLdq(this->m_q1,(q-this->m_q1)/this->m_h )
			//-	this->m_problem.JvdLdv(this->m_q1,(q-this->m_q1)/this->m_h )/this->m_h;
		return true;
	}
};
} // namespace Variational

#endif
