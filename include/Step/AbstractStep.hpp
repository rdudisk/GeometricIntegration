#ifndef DEF_ABSTRACT_STEP
#define DEF_ABSTRACT_STEP

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "include/Common/NOXVector.hpp"

namespace Abstract {

template <typename T_Q, int T_N_EQUATIONS>
class Step
{
public:
	virtual const NOXVector<T_Q::DOF*T_N_EQUATIONS>
	getInitialGuess () = 0;

	virtual bool
	computeF (NOXVector<T_Q::DOF*T_N_EQUATIONS>& f, const NOXVector<T_Q::DOF*T_N_EQUATIONS>& x) = 0;
	
	virtual bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF*T_N_EQUATIONS,T_Q::DOF*T_N_EQUATIONS>& J, const NOXVector<T_Q::DOF*T_N_EQUATIONS>& x) = 0;
};

} // namespace Step

#endif
