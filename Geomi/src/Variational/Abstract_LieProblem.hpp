#ifndef DEF_VARIATIONAL_ABSTRACT_LIEPROBLEM
#define DEF_VARIATIONAL_ABSTRACT_LIEPROBLEM

#include <Eigen/Core>

namespace Variational {
namespace Abstract {

template <typename T_M, typename T_Q, typename T_ALGEBRA>
class LieProblem : public DiscSyst<T_M,T_Q>
{
public:
	LieProblem<T_M,T_Q,T_ALGEBRA> ()
	{ }

	virtual
	~LieProblem<T_M,T_Q,T_ALGEBRA> ()
	{ }
	
	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdv (const T_ALGEBRA) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	JvdLdv (const T_ALGEBRA) = 0;
};

} // namespace Abstract
} // namespace Variational

#endif
