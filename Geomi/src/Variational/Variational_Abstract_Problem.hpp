#ifndef DEF_VARIATIONAL_ABSTRACT_PROBLEM
#define DEF_VARIATIONAL_ABSTRACT_PROBLEM

#include <Eigen/Core>

namespace Variational {
namespace Abstract {

template <typename T_M, typename T_Q>
class Problem : public DiscSyst<T_M,T_Q>
{
public:
	Problem<T_M,T_Q> ()
	{ }

	virtual
	~Problem<T_M,T_Q> ()
	{ }
	
	/**
	 * Compute the derivative of the Lagrangian \f$L\f$ with respect to the position variable \f$q\f$, that is
	 * \f[ \frac{\partial L}{\partial q}(q,v) = \begin{bmatrix}\frac{\partial L}{\partial q_1}\\\vdots\\\frac{\partial L}{\partial q_N}\end{bmatrix}(q,v) \f]
	 */
	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdq (const T_Q, const T_Q) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdv (const T_Q, const T_Q) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	JqdLdq (const T_Q, const T_Q) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	JvdLdq (const T_Q, const T_Q) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	JqdLdv (const T_Q, const T_Q) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	JvdLdv (const T_Q, const T_Q) = 0;
};

} // namespace Abstract
} // namespace Variational

#endif
