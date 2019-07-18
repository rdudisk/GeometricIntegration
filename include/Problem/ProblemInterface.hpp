#ifndef DEF_ABSTRACT_PROBLEM
#define DEF_ABSTRACT_PROBLEM

#include "include/Problem/DiscSystem.hpp"

namespace Abstract {

template <typename T_M, typename T_Q, typename T_LIE_ALGEBRA>
class Problem : public DiscSyst<T_M,T_Q>
{
public:
	Problem<T_M,T_Q,T_LIE_ALGEBRA> ()
	{ }

	virtual
	~Problem<T_M,T_Q,T_LIE_ALGEBRA> ()
	{ }

	virtual bool
	computeA (T_LIE_ALGEBRA& A, const T_Q& x) = 0;

	virtual bool
	computeJacobianA (std::vector<T_LIE_ALGEBRA>& JA, const T_Q& x) = 0;
	
};

} // namespace Abstract

#endif
