#ifndef DEF_SO3_UNCONSTRAINED
#define DEF_SO3_UNCONSTRAINED

#include <cmath>

#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/ellint_1.hpp>

#include "DiscLagSyst.hpp"
#include "SO3Algebra.hpp"
#include "SO3Group.hpp"

#define sgn(x) ((x>=0.0)?1:-1)
#define minus_1_pow(n) ((n%2==1)?-1.0:1.0)

template<typename T>
struct Params {
	Eigen::Matrix<T,3,3> I;
};

using namespace Lie;
using namespace SO3;

template<typename T, typename M>
class DLS_temp : public DiscLagSyst<M,Group<T>,Algebra<T>,Params<T>> {
private:
	int m_exact_initialized = 0;
	int m_orderflag;
	T m_Lnorm;
	T m_w1m;
	T m_w2m;
	T m_w3m;
	T m_I1;
	T m_I2;
	T m_I3;
	T m_wp;
	T m_k;
	T m_eps;
	T m_K;
	T* m_cr;
	T* m_ci;
	T m_A1;
	T m_A2;
	int m_NT;
	Eigen::Matrix<T,3,3> m_P;
	Eigen::Matrix<T,3,3> m_TU0;

	void initialize () {
		
		// First we need to make shure that I2 is the middle one, i.e. I1<=I2<=I3 or I1>=I2>=I3
		
		T	Ix = this->m_params.I(0,0),
			Iy = this->m_params.I(1,1),
			Iz = this->m_params.I(2,2);
		m_P = Eigen::Matrix<T,3,3>::Identity();
		if ((Iy>Ix && Iy>Iz) || (Iy<Ix && Iy<Iz)) {
			m_P << 0, 0,  1, 1, 0, 0, 0, 1, 0;
		}
		Eigen::Matrix<T,3,1> Ivec = m_P*Eigen::Matrix<T,3,1>(Ix,Iy,Iz);
		m_I1 = Ivec[0];
		m_I2 = Ivec[1];
		m_I3 = Ivec[2];
		//std::cout << "I : " << m_I1 << " " << m_I2 << " " << m_I3 << std::endl;

		// Initialisation of the speed limit condition
		// TODO need to work on a better init
		
		Eigen::Matrix<T,3,1> w0 = m_P*(this->m_node[0]).vel().toVector();
		T	w10 = w0[0],
			w20 = w0[1], 
			w30 = w0[2];
		// Override of w0 for testing purpose
		w10 = T(1);
		w20 = T(5);
		w30 = T(6);

		// We need to check the Jacobi condition		

		T	L2 = m_I1*m_I1 * w10*w10  +  m_I2*m_I2 * w20*w20  +  m_I3*m_I3 * w30*w30;
		m_Lnorm = sqrt(L2);
		T	EE = m_I1 * w10*w10  +  m_I2 * w20*w20  +  m_I3 * w30*w30;
		T	Lperp;
		if ((EE>L2/m_I2 && m_I1<m_I3) || (EE<L2/m_I2 && m_I1>m_I3)) {
			m_orderflag = 1;
			T tmp = m_I1;
			m_I1 = m_I3;
			m_I3 = tmp;
			tmp = w10;
			w10 = w30;
			w20 = -w20;
			w30 = tmp;
			Lperp = sqrt(m_I1*m_I1*w10*w10+m_I2*m_I2*w20*w20);
			m_TU0 <<	-Lperp/m_Lnorm,		-m_I2*m_I3*w20*w30/(m_Lnorm*Lperp),	m_I1*m_I3*w10*w30/(m_Lnorm*Lperp),
						0,					-m_I1*w10/Lperp,					-m_I2*w20/Lperp,
						m_I3*w30/m_Lnorm,	-m_I2*w20/m_Lnorm,					m_I1*w10/m_Lnorm;
		} else {
			m_orderflag = 0;
			Lperp = sqrt(m_I1*m_I1*w10*w10+m_I2*m_I2*w20*w20);
			m_TU0 <<	m_I1*m_I3*w10*w30/(m_Lnorm*Lperp),	-m_I2*m_I3*w20*w30/(m_Lnorm*Lperp),	-Lperp/m_Lnorm,
						-m_I2*w20/Lperp,					-m_I1*w10/Lperp,					0,
						m_I1*w10/m_Lnorm,					-m_I2*w20/m_Lnorm,					m_I3*w30/m_Lnorm;
		}
		// Note : since we do not use the initial attitude, B matrix in van Zon article is simply TU0

		// We compute the various coeffs

		m_w1m	=	sgn(w10) * sqrt( (L2-EE*m_I3) / (m_I1*(m_I1-m_I3)) );
		m_w2m	=	-sgn(w20) * sqrt( (L2-EE*m_I3) / (m_I2*(m_I2-m_I3)) );
		m_w3m	=	sgn(w30) * sqrt( (L2-EE*m_I1) / (m_I3*(m_I3-m_I1)) );
		m_wp	=	sgn(m_I2-m_I3) * sgn(w30) * sqrt( (L2-EE*m_I1) * (m_I3-m_I2) / (m_I1*m_I2*m_I3) );
		m_k		=	sqrt( (L2-EE*m_I3) * (m_I1-m_I2)  /  ( (L2-EE*m_I1) * (m_I3-m_I2) )  );
		m_eps	=	boost::math::ellint_1<T>( m_k , asin(w20/m_w2m) );
		m_K		=	boost::math::ellint_1<T>( m_k );
		T CK	=	boost::math::ellint_1<T>( sqrt(1.0-m_k*m_k) );
		T q		=	exp( -M_PI*CK/m_K );
		T eta	=	sgn(w30) * CK - boost::math::ellint_1<T>( sqrt(1.0-m_k*m_k) , asin(m_I3*m_w3m/m_Lnorm) );
		T chi	=	exp(M_PI*eta/m_K);

		T m_A2	=	m_Lnorm/m_I1 + M_PI*m_wp*(chi+1.0)/(2.0*m_K*(chi-1.0));
		int n = 1;
		T dA2;
		do {
			dA2		=	-(M_PI*m_wp/m_K) * (pow(q,2*n)/(1-pow(q,2*n))) * (pow(chi,n)-pow(chi,-n));
			m_A2	+=	dA2;
			n++;
		} while ((dA2 > m_A2*std::numeric_limits<T>::epsilon()) || (dA2 >= std::numeric_limits<T>::min()));
		// pour la condition, voir l'exemple sur https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon

		int m_NT = int(ceil( log(std::numeric_limits<T>::epsilon()*q) / log(q) )) + 1; // +1 compared to Van Zon
		// pareil pour celle-la
		T	i0	=	0,
			r0	=	0;
		m_cr	=	new T[m_NT];
		m_ci	=	new T[m_NT];
		for ( n = 0; n<m_NT; n++ ) {
			m_cr[n]	=   minus_1_pow(n)   * 2 * pow(q,n*(n+1)+0.25) * cosh( ((2*n+1)*M_PI*eta) / (2*m_K) );
			m_ci[n] =	minus_1_pow(n+1) * 2 * pow(q,n*(n+1)+0.25) * sinh( ((2*n+1)*M_PI*eta) / (2*m_K) );
			r0	+=	m_cr[n] * sin( (2*n+1)*M_PI*m_eps / (2*m_K) );
			i0	+=	m_ci[n] * cos( (2*n+1)*M_PI*m_eps / (2*m_K) );
		}
		m_A1	=	atan( i0/r0 ) + ((r0>0)?0:sgn(i0))*M_PI;
		m_exact_initialized = 1;
	}

public:
	virtual void step (const JetSpace<M,Group<T>,Algebra<T>>& js0, JetSpace<M,Group<T>,Algebra<T>>& js1) {
		M h = js1.base()-js0.base();
		js1.pos(js0.pos());
		js1.vel(js0.vel());
	}
	
	JetSpace<M,Group<T>,Algebra<T>> exact_solution (M const& t) {
		if (!m_exact_initialized)
			initialize();

		T sn, cn, dn;
		T wpteps = m_wp*t*m_eps;
		sn = boost::math::jacobi_elliptic<T>(m_k,wpteps,&cn,&dn);
		T	w1 = m_w1m*cn,
			w2 = m_w2m*sn,
			w3 = m_w3m*dn;
		T	Ret1 = 0,
			Imt1 = 0;
		for (int i=0;i<m_NT;i++) {
			Ret1 += m_cr[i] * sin( (2*i+1)*M_PI*wpteps /(2*m_K) );
			Imt1 += m_ci[i] * cos( (2*i+1)*M_PI*wpteps /(2*m_K) );
		}
		T	C = cos( m_A1+m_A2*t ),
			S = sin( m_A1+m_A2*t );
		T cphi = (C*Ret1+S*Imt1) /sqrt(Ret1*Ret1+Imt1*Imt1);
		T sphi = (S*Ret1-C*Imt1) /sqrt(Ret1*Ret1+Imt1*Imt1);
		T Lperp = sqrt( m_I1*m_I1*w1*w1 + m_I2*m_I2*w2*w2 );
		Eigen::Matrix<T,3,3> UT;
		UT <<	m_I1*m_I3*w1*w3/(m_Lnorm*Lperp),	-m_I2*w2/Lperp,	m_I1*w1/m_Lnorm,
				m_I2*m_I3*w2*w3/(m_Lnorm*Lperp),	m_I1*w1/Lperp,	m_I2*w2/m_Lnorm,
				-Lperp/m_Lnorm,				T(0),			m_I3*w3/m_Lnorm;
		if (m_orderflag) {
			UT <<	-Lperp/m_Lnorm,				T(0),			m_I3*w3/m_Lnorm,
					-m_I2*m_I3*w2*w3/(m_Lnorm*Lperp),	-m_I1*w1/Lperp,	-m_I2*w2/m_Lnorm,
					m_I1*m_I3*w1*w3/(m_Lnorm*Lperp),	-m_I2*w2/Lperp,	m_I1*w1/m_Lnorm;
			T tmp = w1;
			w1 = w3;
			w3 = tmp;
			w2 = -w2;
		} 
		Eigen::Matrix<T,3,3> V;
		V <<	cphi, sphi, T(0),
				-sphi, cphi, T(0),
				T(0), T(0), T(1);
		Eigen::Matrix<T,3,1> w;
		w << w1, w2, w3;
		JetSpace<M,Group<T>,Algebra<T>> js;
		Group<T> pos;
		pos.q((m_P.transpose()*UT*V*m_TU0).transpose());
		js.pos(pos);
		Algebra<T> vel;
		std::cout << m_P.transpose() << std::endl;
		std::cout << w;
		Eigen::Matrix<T,3,1> mat_vel(m_P.transpose()*w);
		vel.v(mat_vel);
		js.vel(vel);
		return js;
	}

	void evolve_exact () {
		for (int i=1; i<this->m_node.size()-1; i++) {
			this->m_node[i+1] = exact_solution(this->m_node[i+1].base());
		}
	}
};

typedef float T;
typedef DLS_temp<T,T> DLS;

#endif
