#ifndef DEF_BIVARIATIONAL_ABSTRACT_PROBLEM
#define DEF_BIVARIATIONAL_ABSTRACT_PROBLEM

#include <Eigen/Core>

namespace BiVariational {
namespace Abstract {

template <typename T_M, typename T_Q>
class Problem : public DiscSyst<T_M,T_Q>
{
public:
	Problem ()
	{ }

	virtual
	~Problem ()
	{ }
	
	/**
	 * Compute the derivative of the Lagrangian \f$L\f$ with respect to the
	 * position variable \f$q\f$, that is
	 * \f[
	 *	\frac{\partial L}{\partial q}(q,v) =
	 *		\begin{bmatrix}\frac{\partial L}{\partial q_1}
	 *			\\\vdots\\\frac{\partial L}{\partial q_N}\end{bmatrix}(q,v)
	 *	\f]
	 */
	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdq (const T_Q, const T_Q) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdv0 (const T_Q, const T_Q) = 0;
	
	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdv1 (const T_Q, const T_Q) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	JqdLdq (const T_Q, const T_Q) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	Jv0dLdv0 (const T_Q, const T_Q) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	Jv1dLdv1 (const T_Q, const T_Q) = 0;
};

} // namespace Abstract
} // namespace Variational

#endif
