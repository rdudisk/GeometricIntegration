#ifndef DEF_COMMON_UTILS
#define DEF_COMMON_UTILS

#include <limits>
#include <cmath>

const double BERNOULLI_NUMBERS[] = {
		1.0,			/* 0 */
		-1.0/2.0,	
		1.0/6.0,	
		0.0,
		-1.0/30.0,		/* 5 */
		0.0,
		1.0/42.0,	
		0.0,	
		-1.0/30.0,
		0.0,	
		5.0/66.0,		/* 10 */
		0.0,
		-691.0/2730.0,	
		0.0,	
		7.0/6.0,
		0.0,			/* 15 */
		-3617.0/510.0,	
		0.0,
		43867.0/798.0,	
		0.0,	
		-174611.0/330.0	/* 20 */
	};

// Voir https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/float_comparison.html
template <typename T_SCALAR_TYPE>
bool
isZero (T_SCALAR_TYPE d)
{
	return (std::abs(d)<std::numeric_limits<T_SCALAR_TYPE>::epsilon());
}

template bool isZero (double);


#endif
