#ifndef DEF_RANDOM_NAME
#define DEF_RANDOM_NAME

#include "headers/DiscLagSyst.hpp"

${includes}
${typedefs}

typedef JetSpace<M,Q,TQ> JS;
typedef DiscJetSpace<M,Q,TQ> DJS;

${params}

void step (const JetSpace<M,Q,TQ>& js0, JetSpace<M,Q,TQ>& js1) {
${step_impl}
}

typedef DiscLagSyst<M,Q,TQ,Params> DLS;

#endif
