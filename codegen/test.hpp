#ifndef DEF_RANDOM_NAME
#define DEF_RANDOM_NAME

#include "headers/DiscLagSyst.hpp"

#include "headers/Vec.hpp"

typedef Vec<float,2> vec2;

typedef float M;
typedef vec2 Q;
typedef vec2 TQ;

typedef JetSpace<M,Q,TQ> JS;
typedef DiscJetSpace<M,Q,TQ> DJS;

struct Params {
	float m;
};

void step (const JetSpace<M,Q,TQ>& js0, JetSpace<M,Q,TQ>& js1) {
	M h = js1.base()-js0.base();
	js1.pos(h*js0.vel() + js0.pos() + 2);
	js1.vel(h/(this->m_params.m*pow(h*js0.vel() + js0.pos(), 2)) + js0.vel());
}

typedef DiscLagSyst<M,Q,TQ,Params> DLS;

#endif
