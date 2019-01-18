#ifndef DEF_TEST
#define DEF_TEST

#include "headers/DiscLagSyst.hpp"

#include "headers/Vec.hpp"

typedef Vec<float,2> vec2;

typedef	float M;
typedef	vec2 Q;
typedef	vec2 TQ;

typedef JetSpace<M,Q,TQ> JS;
typedef DiscJetSpace<M,Q,TQ> DJS;

struct Params {
	float m;
};

void step (const JetSpace<M,Q,TQ>& js0, JetSpace<M,Q,TQ>& js1) {
	js1.pos(js0.pos());
	js1.vel(js0.vel());
}

typedef DiscLagSyst<M,Q,TQ,Params> DLS;

#endif
