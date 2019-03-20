#ifndef DEF_DISC_LAG_SYST
#define DEF_DISC_LAG_SYST

#include <iostream>
#include <cmath>
#include <vector>

#include "DiscJetSpace.hpp"

template <typename M, typename Q, typename TQ, typename P>
class DiscLagSyst
: public DiscJetSpace<M,Q,TQ>
{
protected:
	P m_params;

public:
	DiscLagSyst<M,Q,TQ,P> ( )
	{ }

	~DiscLagSyst<M,Q,TQ,P> ( )
	{ }

	P
	params ( ) const {
		return m_params;
	}

	void
	params (P _params) {
		m_params = _params;
	}

	Eigen::Matrix<double,Eigen::Dynamic,1>
	m_dLdq (const Q, const Q);

	Eigen::Matrix<double,Eigen::Dynamic,1>
	m_dLdv (const Q, const Q);

	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
	m_JqdLdq (const Q, const Q);

	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
	m_JvdLdq (const Q, const Q);

	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
	m_JqdLdv (const Q, const Q);

	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
	m_JvdLdv (const Q, const Q);

	//virtual void step (const JetSpace<M,Q,TQ>& js0, JetSpace<M,Q,TQ>& js1) = 0;

	/**
	 * Make the system evolve from initial conditions defined by the JetSpace element `this[0]`
	 * by applying the one step method \p step provided.
	 * The one step method \p step must implement the discrete flow of the Lagrangian associated to
	 * the system.
	 */
	void
	evolve ( void (*step) (const JetSpace<M,Q,TQ>&, JetSpace<M,Q,TQ>&) )
	{
		for (int i=0; i<this->m_node.size()-1; i++) {
			(*step)(this->m_node[i],this->m_node[i+1]);
		}
	}
};

#endif
