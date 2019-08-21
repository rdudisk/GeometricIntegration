#ifndef DEF_VARIATIONAL_IMPLICITEULERSTEPINTERNALS
#define DEF_VARIATIONAL_IMPLICITEULERSTEPINTERNALS

#include <vector>

namespace Variational {

template <typename T_M,
		  typename T_Q,
		  typename T_TQ>
class ImplicitEulerStepInternals : public Abstract::StepInternals<T_M,T_Q,T_TQ>, public ::Abstract::NOXStep<T_Q,1>
{
public:
	ImplicitEulerStepInternals<T_M,T_Q,T_TQ> (Abstract::Problem<T_M,T_Q>& problem)
	:	Abstract::StepInternals<T_M,T_Q,T_TQ>(problem)
	{ }

	/*
	 * T_Q
	posFromVel (T_Q q0, T_TQ v0)
	{
		return q0+this->m_h*v0;
	}
	*/

	const NOXVector<T_Q::DOF>&
	getInitialGuess ()
	{
		NOXVector<T_Q::DOF>* ret = new NOXVector<T_Q::DOF>((1.0+1.0/this->m_h)*this->m_q1-(1.0/this->m_h)*this->m_q0);
		return *ret;
	}

	bool
	computeF (NOXVector<T_Q::DOF>& f, const NOXVector<T_Q::DOF>& q)
	{
		f =	this->m_h *		this->m_problem.dLdq(this->m_q1,(this->m_q1-this->m_q0)/this->m_h)
						+	this->m_problem.dLdv(this->m_q1,(this->m_q1-this->m_q0)/this->m_h)
						-	this->m_problem.dLdv(q,(q-this->m_q1)/this->m_h);
		return true;
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>& J, const NOXVector<T_Q::DOF>& q)
	{
		J =		this->m_problem.JqdLdv(q,(q-this->m_q1)/this->m_h )
			-	this->m_problem.JvdLdv(q,(q-this->m_q1)/this->m_h )/this->m_h;
		return true;
	}
};
} // namespace Variational

#endif
