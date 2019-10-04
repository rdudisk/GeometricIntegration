#ifndef DEF_MULTIVARIATIONAL_ABSTRACT_LIEPROBLEM
#define DEF_MULTIVARIATIONAL_ABSTRACT_LIEPROBLEM

#include <Eigen/Core>

namespace MultiVariational {
namespace Abstract {

template <typename T_Q, typename T_ALGEBRA, int T_N_INDEP>
class LieProblem : public MultiVariational::DiscSyst<T_Q,T_N_INDEP>
{
public:
	LieProblem<T_Q,T_ALGEBRA,T_N_INDEP> ()
	{ }

	virtual
	~LieProblem<T_Q,T_ALGEBRA,T_N_INDEP> ()
	{ }

	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdv (int dir, const T_ALGEBRA) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	JvdLdv (int der_dir, int jac_dir, const T_ALGEBRA) = 0;
};

} // namespace Abstract
} // namespace MultiVariational

#endif
