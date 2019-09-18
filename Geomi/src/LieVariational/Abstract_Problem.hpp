#ifndef DEF_LIEVARIATIONAL_ABSTRACT_PROBLEM
#define DEF_LIEVARIATIONAL_ABSTRACT_PROBLEM

#include <Eigen/Core>

namespace LieVariational {
namespace Abstract {

template <typename T_M, typename T_Q, typename T_ALGEBRA>
class Problem : public DiscSyst<T_M,T_Q>
{
public:
	Problem<T_M,T_Q,T_ALGEBRA> ()
	{ }

	virtual
	~Problem<T_M,T_Q,T_ALGEBRA> ()
	{ }
	
	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdv (const T_ALGEBRA) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	JvdLdv (const T_ALGEBRA) = 0;
};

} // namespace Abstract
} // namespace Variational

#endif
