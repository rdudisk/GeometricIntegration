#ifndef DEF_VARIATIONAL_EULERSTEPINTERNALS
#define DEF_VARIATIONAL_EULERSTEPINTERNALS

#include <vector>

namespace Variational {

template <typename T_M,
		  typename T_Q>
class EulerStepInternals : public Abstract::StepInternals<T_M,T_Q,Abstract::Problem<T_M,T_Q>>, public ::Abstract::NOXStep<T_Q,1>
{
protected:
	using Abstract::StepInternals<T_M,T_Q,Abstract::Problem<T_M,T_Q>>::m_h;
	using Abstract::StepInternals<T_M,T_Q,Abstract::Problem<T_M,T_Q>>::m_q0;
	using Abstract::StepInternals<T_M,T_Q,Abstract::Problem<T_M,T_Q>>::m_q1;
	using Abstract::StepInternals<T_M,T_Q,Abstract::Problem<T_M,T_Q>>::m_problem;

public:
	EulerStepInternals<T_M,T_Q> (Abstract::Problem<T_M,T_Q>& problem)
	:	Abstract::StepInternals<T_M,T_Q,Abstract::Problem<T_M,T_Q>>(problem)
	{ }

	const NOXVector<T_Q::DOF>
	getInitialGuess ()
	{
		NOXVector<T_Q::DOF> ret((1.0+1.0/m_h)*m_q1-(1.0/m_h)*m_q0);
		return ret;
	}

	T_Q
	posFromVel (T_M h, T_Q q0, T_Q v0) const
	{ return q0+h*v0; }

	bool
	computeF (NOXVector<T_Q::DOF>& f, const NOXVector<T_Q::DOF>& q)
	{
		f =				m_problem.dLdv(m_q0,(m_q1-m_q0)/m_h)
			+ m_h *		m_problem.dLdq(m_q1,(q-m_q1)/m_h)
					-	m_problem.dLdv(m_q1,(q-m_q1)/m_h);
		return true;
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>& J, const NOXVector<T_Q::DOF>& q)
	{
		J =		m_problem.JvdLdq(m_q1,(q-m_q1)/m_h)
			-	m_problem.JvdLdv(m_q1,(q-m_q1)/m_h)/m_h;
		return true;
	}
};
} // namespace Variational

#endif
