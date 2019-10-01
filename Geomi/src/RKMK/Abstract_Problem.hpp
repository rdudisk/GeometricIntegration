#ifndef DEF_RKMK_ABSTRACT_PROBLEM
#define DEF_RKMK_ABSTRACT_PROBLEM

//#include "../Common/Common_DiscSyst.hpp"

namespace RKMK {
namespace Abstract {

template <typename T_LIE_ALGEBRA, typename T_M = double>
class Problem : public DiscSyst<T_M,NOXVector<T_LIE_ALGEBRA::DOF>>
{
public:
	Problem<T_LIE_ALGEBRA,T_M> ()
	{ }

	virtual
	~Problem<T_LIE_ALGEBRA,T_M> ()
	{ }

	virtual bool
	computeA (T_LIE_ALGEBRA& A, const NOXVector<T_LIE_ALGEBRA::DOF>& x) = 0;

	virtual bool
	computeJacobianA (std::vector<T_LIE_ALGEBRA>& JA, const NOXVector<T_LIE_ALGEBRA::DOF>& x) = 0;
	
};

} // namespace Abstract
} // namespace RKMK

#endif
