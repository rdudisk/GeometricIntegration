#ifndef DEF_VARIATIONAL_MIDPOINTSTEPINTERNALS
#define DEF_VARIATIONAL_MIDPOINTSTEPINTERNALS

#include <vector>

namespace Variational {

template <typename T_M,
		  typename T_Q,
		  typename T_TQ>
class MidpointStepInternals : public Abstract::StepInternals<T_M,T_Q,T_TQ>, public ::Abstract::NOXStep<T_Q,1>
{
public:
	MidpointStepInternals<T_M,T_Q,T_TQ> (Abstract::Problem<T_M,T_Q>& problem)
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
		std::cout << "q0+q1/2 = \n" << (this->m_q0+this->m_q1)/2.0 << std::endl;
		std::cout << "q1-q0/h = \n" << (this->m_q1-this->m_q0)/this->m_h << std::endl;
		std::cout << "q1+q/2 = \n" << (this->m_q1+q)/2.0 << std::endl;
		std::cout << "q-q1/h = \n" << (q-this->m_q1)/this->m_h << std::endl;
		f = (this->m_h/2.0) * 	this->m_problem.dLdq((this->m_q0+this->m_q1)/2.0,(this->m_q1-this->m_q0)/this->m_h)
							+	this->m_problem.dLdv((this->m_q0+this->m_q1)/2.0,(this->m_q1-this->m_q0)/this->m_h)
		  + (this->m_h/2.0) *	this->m_problem.dLdq((this->m_q1+q)/2.0,(q-this->m_q1)/this->m_h)
							-	this->m_problem.dLdv((this->m_q1+q)/2.0,(q-this->m_q1)/this->m_h);
		return true;
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>& J, const NOXVector<T_Q::DOF>& q)
	{
		NOXVector<T_Q::DOF> pos = (q+this->m_q1)/2.0;
		NOXVector<T_Q::DOF> vel = (q-this->m_q1)/this->m_h;
		J = (this->m_h/4.0) *	this->m_problem.JqdLdq(pos,vel)
							+	this->m_problem.JvdLdq(pos,vel)/2.0
							-	this->m_problem.JqdLdv(pos,vel)/2.0
							-	this->m_problem.JvdLdv(pos,vel)/this->m_h;
		return true;
	}
};
} // namespace Variational

#endif
