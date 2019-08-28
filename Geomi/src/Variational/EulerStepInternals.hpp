#ifndef DEF_VARIATIONAL_EULERSTEPINTERNALS
#define DEF_VARIATIONAL_EULERSTEPINTERNALS

#include <vector>

namespace Variational {

template <typename T_M,
		  typename T_Q,
		  typename T_TQ>
class EulerStepInternals : public Abstract::StepInternals<T_M,T_Q,T_TQ>, public ::Abstract::NOXStep<T_Q,1>
{
public:
	EulerStepInternals<T_M,T_Q,T_TQ> (Abstract::Problem<T_M,T_Q>& problem)
	:	Abstract::StepInternals<T_M,T_Q,T_TQ>(problem)
	{ }

	const NOXVector<T_Q::DOF>
	getInitialGuess ()
	{
		NOXVector<T_Q::DOF> ret((1.0+1.0/this->m_h)*this->m_q1-(1.0/this->m_h)*this->m_q0);
		return ret;
	}

	bool
	computeF (NOXVector<T_Q::DOF>& f, const NOXVector<T_Q::DOF>& q)
	{
		f =					this->m_problem.dLdv(this->m_q0,(this->m_q1-this->m_q0)/this->m_h)
			+ this->m_h *	this->m_problem.dLdq(this->m_q1,(q-this->m_q1)/this->m_h)
						-	this->m_problem.dLdv(this->m_q1,(q-this->m_q1)/this->m_h);
		return true;
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>& J, const NOXVector<T_Q::DOF>& q)
	{
		J =		this->m_problem.JvdLdq(this->m_q1,(q-this->m_q1)/this->m_h )
			-	this->m_problem.JvdLdv(this->m_q1,(q-this->m_q1)/this->m_h )/this->m_h;
		return true;
	}
};
} // namespace Variational

#endif
