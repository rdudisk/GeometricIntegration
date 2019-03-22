#ifndef DEF_VAR_INTEGRATOR_DEFS
#define DEF_VAR_INTEGRATOR_DEFS

#include <vector>

#include "DiscLagSyst.hpp"

template<typename M, typename Q, typename TQ, typename P>
struct var_integrator_params
{
	M h;
	std::vector<Q> pos;
	DiscLagSyst<M,Q,TQ,P> *syst;
	void* additional_params;
};

#endif
