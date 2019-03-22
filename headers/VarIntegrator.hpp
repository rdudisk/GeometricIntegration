#ifndef DEF_VAR_INTEGRATOR
#define DEF_VAR_INTEGRATOR

#include <iostream>
#include <cmath>
#include <vector>

#include <typeinfo>

#include "VarIntegratorDefs.hpp"
#include "DiscLagSyst.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "Galerkin.hpp"

int
print_state (size_t iter, gsl_multiroot_fsolver *s)
{
	printf("iter = %3u x = % .3f % .3f "
			"f(x) = % .3e % .3e\n",
			iter,
			gsl_vector_get(s->x,0),
			gsl_vector_get(s->x,1),
			gsl_vector_get(s->f,0),
			gsl_vector_get(s->f,1));
}

int
print_state_fdf (size_t iter, gsl_multiroot_fdfsolver *s)
{
	printf("iter = %3u x = % .3f % .3f "
			"f(x) = % .3e % .3e\n",
			iter,
			gsl_vector_get(s->x,0),
			gsl_vector_get(s->x,1),
			gsl_vector_get(s->f,0),
			gsl_vector_get(s->f,1));
}

template <typename M, typename Q, typename TQ, typename P>
class VarIntegrator
{
private:
	/* WARNING: there seems to be something wrong with gsl_vector_set when the vector representing the DEL is not an Eigen::...<double>
	 * This is also the case when computing m_dLdq and m_dLdv internally with Q of scalar type other than double (at least float tested)
	 * and casting it with Eigen::...::cast<double> beofre the return.
	 * Couldn't find the cause of the problem after hours of search, so for now the recommandation is to use double type for scalar representation
	 * within Q (and M).*/

	/* Euler internal functions for solving DEL
	 */
	static int
	f_del_euler (const gsl_vector* x_vec, void* pp, gsl_vector* f)
	{
		var_integrator_params<M,Q,TQ,P>* p = (struct var_integrator_params<M,Q,TQ,P>*) pp;
		M h = p->h;
		const Q x0 = p->pos[0];
		const Q x1 = p->pos[1];
		DiscLagSyst<M,Q,TQ,P>* syst = p->syst;

		const size_t dof = syst->dof();

		Q x(Q::cast_from_gsl_vector(x_vec));

		Eigen::Matrix<double,Eigen::Dynamic,1> del =
			      syst->m_dLdv( x0, (x1-x0)/h )
			+ h * syst->m_dLdq( x1, (x-x1)/h )
			    - syst->m_dLdv( x1, (x-x1)/h );
		
		for (size_t i=0; i<dof; i++)
			gsl_vector_set(f,i,del(i));

		return GSL_SUCCESS;
	}

	static int
	df_del_euler (const gsl_vector* x_vec, void* pp, gsl_matrix* J)
	{
		var_integrator_params<M,Q,TQ,P>* p = (struct var_integrator_params<M,Q,TQ,P>*) pp;
		M h = p->h;
		const Q x0 = p->pos[0];
		const Q x1 = p->pos[1];
		DiscLagSyst<M,Q,TQ,P>* syst = p->syst;

		const size_t dof = syst->dof();

		Q x(Q::cast_from_gsl_vector(x_vec));

		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Jdel =
			  syst->m_JvdLdq( x1, (x-x1)/h )
			- syst->m_JvdLdv( x1, (x-x1)/h )/h;
	
		size_t i, j;
		for (i=0; i<dof; i++) {
			for (j=0; j<dof; j++)
				gsl_matrix_set(J,i,j,Jdel(i,j));
		}

		return GSL_SUCCESS;
	}

	static int
	fdf_del_euler (const gsl_vector* x_vec, void* pp, gsl_vector* f, gsl_matrix* J)
	{
		var_integrator_params<M,Q,TQ,P>* p = (struct var_integrator_params<M,Q,TQ,P>*) pp;
		M h = p->h;
		const Q x0 = p->pos[0];
		const Q x1 = p->pos[1];
		DiscLagSyst<M,Q,TQ,P>* syst = p->syst;

		const size_t dof = syst->dof();

		Q x(Q::cast_from_gsl_vector(x_vec));

		Eigen::Matrix<double,Eigen::Dynamic,1> del =
			      syst->m_dLdv( x0, (x1-x0)/h )
			+ h * syst->m_dLdq( x1, (x-x1)/h )
			    - syst->m_dLdv( x1, (x-x1)/h );
		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Jdel =
			  syst->m_JvdLdq( x1, (x-x1)/h )
			- syst->m_JvdLdv( x1, (x-x1)/h )/h;
	
		size_t i, j;
		for (i=0; i<dof; i++) {
			for (j=0; j<dof; j++)
				gsl_matrix_set(J,i,j,Jdel(i,j));
			gsl_vector_set(f,i,del(i));
		}
		return GSL_SUCCESS;
	}


	/* Midpoint internal functions for solving DEL */
	static int
	f_del_midpoint (const gsl_vector* x_vec, void* pp, gsl_vector* f)
	{
		var_integrator_params<M,Q,TQ,P>* p = (struct var_integrator_params<M,Q,TQ,P>*) pp;
		M h = p->h;
		const Q x0 = p->pos[0];
		const Q x1 = p->pos[1];
		DiscLagSyst<M,Q,TQ,P>* syst = p->syst;

		const size_t dof = syst->dof();

		Q x(Q::cast_from_gsl_vector(x_vec));

		const Q x01 = (x0+x1)/2.0;
		const Q dx01 = (x1-x0)/h;
		const Q x12 = (x1+x)/2.0;
		const Q dx12 = (x-x1)/h;
		Eigen::Matrix<double,Eigen::Dynamic,1> del =
			  h*syst->m_dLdq( x01, dx01 )/2.0
			+   syst->m_dLdv( x01, dx01 )
			+ h*syst->m_dLdq( x12, dx12 )/2.0
			-	syst->m_dLdv( x12, dx12 );

		for (size_t i=0; i<dof; i++)
			gsl_vector_set(f,i,del(i));
		 
		return GSL_SUCCESS;
	}

	static int
	df_del_midpoint (const gsl_vector* x_vec, void* pp, gsl_matrix* J)
	{
		var_integrator_params<M,Q,TQ,P>* p = (struct var_integrator_params<M,Q,TQ,P>*) pp;
		M h = p->h;
		const Q x0 = p->pos[0];
		const Q x1 = p->pos[1];
		DiscLagSyst<M,Q,TQ,P>* syst = p->syst;

		const size_t dof = syst->dof();

		Q x(Q::cast_from_gsl_vector(x_vec));

		const Q x12 = (x1+x)/2.0;
		const Q dx12 = (x-x1)/h;
		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Jdel =
			h*syst->m_JqdLdq( x12, dx12 )/4.0
			+ syst->m_JvdLdq( x12, dx12 )/2.0
			- syst->m_JqdLdv( x12, dx12 )/2.0
			- syst->m_JvdLdv( x12, dx12 )/h;
	
		size_t i, j;
		for (i=0; i<dof; i++) {
			for (j=0; j<dof; j++)
				gsl_matrix_set(J,i,j,Jdel(i,j));
		}
		return GSL_SUCCESS;
	}

	static int
	fdf_del_midpoint (const gsl_vector* x_vec, void* pp, gsl_vector* f, gsl_matrix* J)
	{
		var_integrator_params<M,Q,TQ,P>* p = (struct var_integrator_params<M,Q,TQ,P>*) pp;
		M h = p->h;
		const Q x0 = p->pos[0];
		const Q x1 = p->pos[1];
		DiscLagSyst<M,Q,TQ,P>* syst = p->syst;

		const size_t dof = syst->dof();

		Q x(Q::cast_from_gsl_vector(x_vec));

		const Q x01 = (x0+x1)/2.0;
		const Q dx01 = (x1-x0)/h;
		const Q x12 = (x1+x)/2.0;
		const Q dx12 = (x-x1)/h;
		Eigen::Matrix<double,Eigen::Dynamic,1> del =
			  h*syst->m_dLdq( x01, dx01 )/2.0
			+   syst->m_dLdv( x01, dx01 )
			+ h*syst->m_dLdq( x12, dx12 )/2.0
			-   syst->m_dLdv( x12, dx12 );
		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Jdel =
			h*syst->m_JqdLdq( x12, dx12 )/4.0
			+ syst->m_JvdLdq( x12, dx12 )/2.0
			- syst->m_JqdLdv( x12, dx12 )/2.0
			- syst->m_JvdLdv( x12, dx12 )/h;

		size_t i, j;
		for (i=0; i<dof; i++) {
			for (j=0; j<dof; j++)
				gsl_matrix_set(J,i,j,Jdel(i,j));
			gsl_vector_set(f,i,del(i));
		}
		return GSL_SUCCESS;
	}

public:
	DiscLagSyst<M,Q,TQ,P> *m_syst;

	void
	onestep_method_f (  int (*f_del) (const gsl_vector*, void*, gsl_vector*)  )
	{
		const size_t dof = this->m_syst->dof();

		const gsl_multiroot_fsolver_type* TT;
		gsl_multiroot_fsolver* s;

		TT = gsl_multiroot_fsolver_hybrids;
		s = gsl_multiroot_fsolver_alloc(TT,dof);

		struct var_integrator_params<M,Q,TQ,P> pp;

		Q q0,q1;
		M h;

		gsl_multiroot_function f;
		gsl_vector* x = gsl_vector_alloc(dof);

		int status;
		size_t iter;

		for (int i=1; i<m_syst->size()-1; i++) {
			iter = 0;

			q0  = this->m_syst->pos(i-1);
			q1  = this->m_syst->pos(i);
			h   = this->m_syst->base(i) - this->m_syst->base(i-1);

			std::vector<Q> pos;
			pos.push_back (q0);
			pos.push_back (q1);
			pp  = { h, pos, this->m_syst };
			f   = { f_del , dof , &pp };

			gsl_vector_set(x,0,q1(0));
			gsl_vector_set(x,1,q1(1));

			gsl_multiroot_fsolver_set(s,&f,x);

			std::cout << "Step " << i << ":\n";

			do {
				iter++;
				status = gsl_multiroot_fsolver_iterate(s);
				print_state(iter,s);
				if (status) break;
				status = gsl_multiroot_test_residual(s->f,1e-7);
			} while (status == GSL_CONTINUE  &&  iter < 50);

			std::cout << "Status : " << gsl_strerror(status) << "\n";
			std::cout << "Number of iterations : " << iter << "\n";
			std::cout << "\n";

			this->m_syst->pos(i+1,Q::cast_from_gsl_vector(s->x));
		}

		gsl_multiroot_fsolver_free(s);
		gsl_vector_free(x);
	}

	void
	onestep_method_fdf ( int (*f_del)	(const gsl_vector*, void*, gsl_vector*),
						 int (*df_del)	(const gsl_vector*, void*, gsl_matrix*),
						 int (*fdf_del)	(const gsl_vector*, void*, gsl_vector*, gsl_matrix*) )
	{
		const size_t dof = this->m_syst->dof();

		const gsl_multiroot_fdfsolver_type* TT;
		gsl_multiroot_fdfsolver* s;

		TT = gsl_multiroot_fdfsolver_gnewton;
		s = gsl_multiroot_fdfsolver_alloc(TT,dof);

		struct var_integrator_params<M,Q,TQ,P> pp;

		Q q0,q1;
		M h;

		gsl_multiroot_function_fdf f;
		gsl_vector* x = gsl_vector_alloc(dof);

		int status;
		size_t iter;

		for (int i=1; i<m_syst->size()-1; i++) {
			iter = 0;

			q0  = this->m_syst->pos(i-1);
			q1  = this->m_syst->pos(i);
			h   = this->m_syst->base(i)-this->m_syst->base(i-1);

			std::vector<Q> pos;
			pos.push_back (q0);
			pos.push_back (q1);
			pp  = { h, pos, this->m_syst };
			f   = { f_del, df_del, fdf_del, dof , &pp };

			gsl_vector_set(x,0,q1(0));
			gsl_vector_set(x,1,q1(1));

			gsl_multiroot_fdfsolver_set(s,&f,x);

			std::cout << "Step " << i << ":\n";

			do {
				iter++;
				status = gsl_multiroot_fdfsolver_iterate(s);
				print_state_fdf(iter,s);
				if(status) break;
				status = gsl_multiroot_test_residual(s->f,1e-7);
			} while (status == GSL_CONTINUE  &&  iter < 50);

			std::cout << "Status : " << gsl_strerror(status) << "\n";
			std::cout << "Number of iterations : " << iter << "\n";
			std::cout << "\n";

			this->m_syst->pos(i+1,Q::cast_from_gsl_vector(s->x));
		}

		gsl_multiroot_fdfsolver_free(s);
		gsl_vector_free(x);
	}

	void
	midpoint ( )
	{
		onestep_method_f (this->f_del_midpoint);
	}

	void
	midpoint_fdf ( )
	{
		onestep_method_fdf (this->f_del_midpoint, this->df_del_midpoint, this->fdf_del_midpoint);
	}

	void
	euler ( )
	{
		onestep_method_f (this->f_del_euler);
	}

	void
	euler_fdf ( )
	{
		onestep_method_fdf (this->f_del_euler, this->df_del_euler, this->fdf_del_euler);
	}
};

#endif
