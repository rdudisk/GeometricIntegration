#ifndef DEF_MULTIVARIATIONAL_ABSTRACT_LIEPROBLEM
#define DEF_MULTIVARIATIONAL_ABSTRACT_LIEPROBLEM

#include <Eigen/Core>

namespace MultiVariational {
namespace Abstract {

template <typename T_M, typename T_Q, typename T_ALGEBRA, int T_N_INDEP>
class LieProblem : public DiscSyst<T_M,T_Q>
{
public:
	LieProblem<T_M,T_Q,T_ALGEBRA> ()
	{ }

	virtual
	~LieProblem<T_M,T_Q,T_ALGEBRA> ()
	{ }
	
	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdv (int direction, const T_ALGEBRA) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	JvdLdv (int derivative_direction, int jacobian_direction, const T_ALGEBRA) = 0;
};

} // namespace Abstract
} // namespace MultiVariational

#endif
