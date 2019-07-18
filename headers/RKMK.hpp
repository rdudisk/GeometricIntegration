#ifndef DEF_RKMK
#define DEF_RKMK

#include <iostream>
#include <cmath>
#include <vector>

#include "DiscJetSpace.hpp"

template <typename Q, typename TQ, typename GRP, typename ALG>
class RKMK
{
private:

public:
	DiscJetSpace<M,Q,TQ> *m_space;

	// On suppose que la methode est explicite !!
	void 
	rkmk (const size_t s, const float* a, const float* b)
	{
		Q q0;
		float h;

		size_t n,i,j;
		ALG aksum, bksum;
		std::vector<ALG> k(s);

		size_t sum_order;
		sum_order = max(1,s-2);

		for (n=1; n<m_syst->size()-1; n++) {

			q0  = this->m_syst->pos(n-1);
			h   = this->m_syst->base(n)-this->m_syst->base(n-1);

			for (i=0; i<s; i++) {
				aksum = ALG::Zero();
				for (j=0; j<s; j++) {
					aksum += a[i*s+j]*k[j];
				}
				k[i] = f(h*aksum);
			}

			bksum = ALG::Zero();
			for (i=0; i<s; i++) {
				bksum += b[i]*k[i];
			}
			ALG tildomega = h*bksum;

			this->m_syst->pos(n,tildomega.exp()*q0);
		}
	}

	void
	rkmk4 ( )
	{
		const float[] a = { \
			0.0,0.0,0.0,0.0, \
			0.5,0.0,0.0,0.0, \
			0.0,0.5,0.0,0.0, \
			0.0,0.0,1.0,0.0};
		const float[] b = {1.0/6.0,2.0/6.0,2.0/6.0,1.0/6.0};
		this->rkmk(4,a,b);
	}
};

#endif
