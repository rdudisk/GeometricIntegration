#ifndef DEF_BIVARIATIONAL_ABSTRACT_LIEPROBLEM
#define DEF_BIVARIATIONAL_ABSTRACT_LIEPROBLEM

#include <Eigen/Core>

namespace BiVariational {
namespace Abstract {

template <typename T_SCALAR, typename T_Q, typename T_VEL>
class LieProblem : public DiscSyst<T_SCALAR,T_Q,T_VEL>
{
public:
	LieProblem ()
	{ }

	virtual
	~LieProblem ()
	{ }
	
	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdv0 (const T_VEL) = 0;

	virtual Eigen::Matrix<double,T_Q::DOF,1>
	dLdv1 (const T_VEL) = 0;

	//virtual Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>
	//JvdLdv (const T_VEL) = 0;
};

} // namespace Abstract
} // namespace Variational

#endif
