#ifndef DEF_TEST
#define DEF_TEST

#include "headers/DiscLagSyst.hpp"
#include "headers/VarIntegrator.hpp"

#include "headers/Vec.hpp"

typedef Vec<double,2> vec2;

typedef	double	M;
typedef	vec2	Q;
typedef	vec2	TQ;

typedef JetSpace<M,Q,TQ>		JS;
typedef DiscJetSpace<M,Q,TQ>	DJS;

struct Params
{
	float m;
};

void
step (const JetSpace<M,Q,TQ>& js0, JetSpace<M,Q,TQ>& js1)
{
	js1.pos(js0.pos());
	js1.vel(js0.vel());
}

typedef DiscLagSyst<M,Q,TQ,Params>		DLS;
typedef VarIntegrator<M,Q,TQ,Params>	VI;

template<>
Eigen::Matrix<double,Eigen::Dynamic,1>
DLS::m_dLdq (const Q q0, const Q q1)
{
	Eigen::Matrix<double,Eigen::Dynamic,1> res(Q::dof());
	res = -q0.cwiseInverse()/2.0;
	return res;
}

template<>
Eigen::Matrix<double,Eigen::Dynamic,1>
DLS::m_dLdv (const Q q0, const Q q1)
{
	Eigen::Matrix<double,Eigen::Dynamic,1> res(Q::dof());
	res = q1;
	return res;
}

template<>
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
DLS::m_JqdLdq (const Q q0, const Q q1)
{
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> res = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(Q::dof(),Q::dof());
	for (int i=0; i<Q::dof(); i++) {
		res(i,i) = 1.0/(4.0*q0(i)*q0(i));
	}
	return res;
}

template<>
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
DLS::m_JvdLdq (const Q q0, const Q q1)
{
	return Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(Q::dof(),Q::dof());
}

template<>
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
DLS::m_JqdLdv (const Q q0, const Q q1)
{
	return Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Zero(Q::dof(),Q::dof());
}

template<>
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
DLS::m_JvdLdv (const Q q0, const Q q1)
{
	return Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>::Identity(Q::dof(),Q::dof());
}

#endif
