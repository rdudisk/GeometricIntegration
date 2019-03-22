#ifndef DEF_GALERKIN
#define DEF_GALERKIN

#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_vector.h>

#include "../libs/fastgl/fastgl.h"

#include "VarIntegratorDefs.hpp"

struct galerkin_params
{
	size_t poly_deg;
	size_t gauss_deg;
	std::vector<double> ctrl;
};

namespace fastgl
{
	/* Structures storing quadrature points in x space rather than the original theta space.
	 */
	struct NodeWeight {
		double c, w;
		// c in [-1,1]

		NodeWeight ( )
		{ }

		NodeWeight (double _c, double _w)
		: c (_w), w (_w)
		{ }
	};
} // namepace fastgl

/*
 * Workaround to have static member in class Galerkin always initialised without
 * the need of explicitely calling a static init function.
 */
class GaussLegendrePairs
{
private:
	/* Stores the values of nodes and weights for a Gauss-Legendre quadrature.
	 * Pairs corresponding to an order r quadrature are stored in the vector m_gl_pairs[deg-1],
	 * and are stored in increasing order, from index 0 to deg-1.
	 * They are evaluated for a quadrature in [-1,1].
	 */
	static
	std::vector<std::vector<fastgl::NodeWeight>> m_gl_pairs;

public:
	/*
	 * Initialize nodes and weights for Gauss-Legendre quadrature up to order 10.
	 */
	void
	init ( )
	{
		const size_t max = 10;
		GaussLegendrePairs::m_gl_pairs = std::vector<std::vector<fastgl::NodeWeight>>(max);
		size_t deg, j;
		for (deg=0; deg<max; deg++) {
			GaussLegendrePairs::m_gl_pairs[deg] = std::vector<fastgl::NodeWeight>(deg);
			for (j=0; j<=deg; j++) {
				fastgl::QuadPair q = fastgl::GLPair (deg+1, j+1); // in theta space
				GaussLegendrePairs::m_gl_pairs[deg][j] = fastgl::NodeWeight (q.x(), q.weight); // in x space
			}
		}
	}

	GaussLegendrePairs ( )
	{
		GaussLegendrePairs::init();
	}

	/*
	 * deg >= 0, k in [1..deg]
	 */
	static fastgl::NodeWeight
	gl_pair (size_t deg, size_t k)
	{
		return GaussLegendrePairs::m_gl_pairs[deg-1][k-1];
	}

	static std::vector<fastgl::NodeWeight>
	gl_pair_list (size_t deg) 
	{
		return GaussLegendrePairs::m_gl_pairs[deg-1];
	}
};

template <typename M, typename Q, typename TQ, typename P>
class Galerkin
{
private:
	static GaussLegendrePairs m_gl;

public:
	/*
	 * Computes all Lagrange polynomials of degree deg evaluated at s, as well as their derivatives.
	 * They are stored in a vector of size 2*deg, first all polynomes evaluations then all derivatives evaluations.
	 */
	std::vector<double>
	p_dp_lag (const size_t deg, const std::vector<double> ctrl, const double s)
	{
		std::vector<double> p_lag (2*deg, 0.0); // evalation des polynomes et des derivees

		// pour la definition de prod, voir cahier 2 p. 84
		std::vector<double> prod (deg*deg, 1.0); // 'deg' espaces non utilises mais simplifie l'acces par indices

		size_t k, nu, i;
		double c;
		for (nu=0; nu<deg; nu++) {
			for (i=0; i<deg; i++) {
				if (i!=nu) {
					c = (s-ctrl[i]) / (ctrl[nu]-ctrl[i]);
					for (k=0; k<deg; k++) {
						if (k!=nu && k!=i) {
							prod[k*deg+nu] *= c;
						}
					}
				}
			}
		}

		for (nu=0; nu<deg; nu++) {
			p_lag[nu] = prod[nu] * (s-ctrl[0]) / (ctrl[nu]-ctrl[0]); // prod[nu] equiv prod(0,nu)
			for (k=0; k<deg; k++) {
				if (k!=nu) {
					p_lag[deg+nu] += prod[k*deg+nu] / (ctrl[nu]-ctrl[k]);
				}
			}
		}

		return p_lag;
	}

	/* Compute the discrete Lagrangian for a galerkin interpolation with Lagrange polynomials and Gauss-Legendre quadrature.
	 */
	double
	disc_lag_galerkin (	double (*lag) (const Q, const TQ),
						const std::vector<Q> q_coeffs,
						const M h,
						const std::vector<double> control_points_p_lag,
						const size_t poly_deg,
						const size_t gauss_deg)
	{
		double sum = 0;

		size_t i, j;
		double c;
		Q q;
		TQ v;
		fastgl::NodeWeight quad;
		for (i=0; i<gauss_deg; i++) {
			q = Q::Zero ( );
			v = TQ::Zero ( );
			quad = Galerkin::m_gl.gl_pair (gauss_deg, i+1);
			c = (1.0+quad.c) / 2.0; // c now takes value in [0,1]
			std::vector<double> p_dp_val = p_dp_lag (poly_deg, control_points_p_lag, c); // all polynomials and derivatives evaluated at c
			for (j=0; j<poly_deg; j++) {
				q += q_coeffs[j] * p_dp_val[j];
				v += q_coeffs[j] * p_dp_val[poly_deg+j];
			}
			sum += quad.w * lag (q, v/h);
		}
		sum *= h/2.0; // multiplication by h/2 to impact the weight with the change of interval from [-1,1] to [0,h]

		return sum;
	}

	static int
	f_del_galerkin (const gsl_vector* x_vec, void* pp, gsl_vector* f)
	{
		/*
		 * convention : s = poly_deg
		 *
		 * curr_pos :	std::vector<Q> de taille poly_deg+1
		 *				current position
		 *					contient des valeurs des points de controle q_0,...,q_s
		 * former_pos : std::vector<Q> de taille poly_deg+1
		 *				former position
		 *					meme chose que curr_pos pour le pas de temps precedent
		 * del :		Eigen::Matrix<double,dof,poly_deg>
		 *				discrete euler lagrange
		 *					equations DEL pour j in [1,poly_deg]
		 *					col 1 : D_{s+1}(former) + D_1(current)
		 *					col j : D_j(current)
		 * p_dp_val :	std::vector<std::vector<double>> size (gauss_deg by (poly_deg+1))
		 *				polynome and derivative values
		 *					p_dp_val [i][j] = l_{j,s}(c_{h,i}/h)
		 * dlag_val :	stc::vector<Eigen::Matrix<double,dof,1>> size (2*gauss_deg)
		 *				derivatives of lagrangian
		 *					dlag_val [2*i]		= dLdq (q_d(c_{h,i}), v_d(c_{h,i}))
		 *					dlag_val [2*i+1]	= dLdv (q_d(c_{h,i}), v_d(c_{h,i}))
		 */
		var_integrator_params<M,Q,TQ,P>* p = (struct var_integrator_params<M,Q,TQ,P>*) pp;
		const M h = p->h;
		const std::vector<Q> former_pos = p->pos;
		DiscLagSyst<M,Q,TQ,P>* syst = p->syst;
		galerkin_params* gp = (struct galerkin_params*) (p->additional_params);

		const size_t dof = syst->dof;
		const size_t poly_deg = gp->poly_deg;
		const size_t gauss_deg = gp->gauss_deg;
		const std::vector<double> ctrl = gp->ctrl;

		std::vector<Q> curr_pos(poly_deg+1);
		gsl_vector* tmp_vec = gsl_vector_alloc(dof);
		size_t i, j, k, nu;

		curr_pos[0] = former_pos[poly_deg];
		for (i=1; i<=poly_deg; i++) {
			for (j=0; j<dof; j++) {
				gsl_vector_set (tmp_vec, j, gsl_vector_get (x_vec, (i-1)*dof+j));
			}
			curr_pos[i] = Q::cast_from_gsl_vector (x_vec);
		}

		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> del(dof,poly_deg);
		
		std::vector<Eigen::Matrix<double,Eigen::Dynamic,1>> dlag_val (2*gauss_deg);
		std::vector<std::vector<double>> p_dp_val (gauss_deg);
		double c;
		std::vector<fastgl::NodeWeight> quad = Galerkin::m_gl.gl_pair_list (gauss_deg);
		Q q;
		TQ v;

		for (i=0; i<gauss_deg; i++) {
			c = (1.0+quad[i].c) / 2.0;
			p_dp_val[i] = p_dp_lag (poly_deg, ctrl, c);
			q = Q::Zero();
			v = TQ::Zero();
			for (nu=0; nu<=poly_deg; nu++) {
				q += curr_pos[nu] * p_dp_val[i][nu];
				v += curr_pos[nu] * p_dp_val[i][poly_deg+nu];
			}
			v /= h;
			dlag_val[2*i]	= syst->dLdq (q, v);
			dlag_val[2*i+1] = syst->dLdv (q, v);
		}
		// Calcul des DEL [0..s-1]
		for (nu=0; nu<poly_deg; nu++) {
			del.col(j) = Eigen::Matrix<double,Eigen::Dynamic,1>::Zero(dof);
			for (i=0; i<gauss_deg; i++) {
				del.col(j) += (h*quad[i].w/2.0) * (dlag_val[2*i]*p_dp_val[i][nu] + dlag_val[2*i+1]*p_dp_val[i][poly_deg+nu]/h);
			}
		}
		// Calcul de DEL_s
		q = Q::Zero();
		v = TQ::Zero();
		for (nu=0; nu<=poly_deg; nu++) {
			q += former_pos[nu] * p_dp_val[poly_deg][nu];
			v += former_pos[nu] * p_dp_val[poly_deg][poly_deg+nu];
		}
		v /= h;
		for (i=0; i<gauss_deg; i++) {
			del.col(0) += (h*quad[i].w/2.0) * (syst->dLdq (q, v)*p_dp_val[i][poly_deg] + syst->dLdv (q, v)*p_dp_val[i][2*poly_deg]/h);
		}

		for (i=0; i<dof; i++) {
			for (nu=0; nu<poly_deg; nu++) {
				gsl_vector_set (f, i+nu*dof, del (i,nu));
			}
		}

		return GSL_SUCCESS;
	}
};

#endif
