#ifndef DEF_RKMK_STEPINTERNALS
#define DEF_RKMK_STEPINTERNALS

#include <vector>

namespace RKMK {

template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES>
class StepInternals
{
protected:
	T_Q m_y0;
	T_M m_h;

	T_LIE_ALGEBRA m_solution;

	Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& m_interface;

	double			m_a_coeffs	[T_N_INTERNAL_STAGES * T_N_INTERNAL_STAGES];
	double			m_b_coeffs	[T_N_INTERNAL_STAGES];
	T_LIE_ALGEBRA	m_k			[T_N_INTERNAL_STAGES];
	unsigned int	m_order_q;

public:
	StepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface)
	:	m_interface(interface)
	{ setTruncatureOrder(); }

	static StepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>&
	newFromProblem (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface)
	{
		StepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>* res = new StepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>(interface);
		return *res;
	}

	void
	operator= (const StepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>& other)
	{
		// TODO: Cod'e avec le cul, il faut verifier tout ca.
		m_y0		= other.m_y0;
		m_h			= other.m_h;
		m_solution	= other.m_solution;
		m_interface = other.m_interface;
		std::copy(std::begin(other.m_a_coeffs),std::end(other.m_a_coeffs),std::begin(m_a_coeffs));
		std::copy(std::begin(other.m_b_coeffs),std::end(other.m_b_coeffs),std::begin(m_b_coeffs));
		std::copy(std::begin(other.m_k),std::end(other.m_k),std::begin(m_k));
		m_order_q	= other.m_order_q;
	}

	double
	a_coeffs (int i, int j) const
	{ return m_a_coeffs[i*T_N_INTERNAL_STAGES+j]; }

	void
	setData (T_M h_var, T_Q y0_var)
	{
		m_h = h_var;
		m_y0 = y0_var;
	}

	bool
	setCoeffs (std::vector<double> a, std::vector<double> b)
	{
		if ((a.size()!=T_N_INTERNAL_STAGES*T_N_INTERNAL_STAGES) || (b.size()!=T_N_INTERNAL_STAGES))
			return false;
		for (int i=0; i<T_N_INTERNAL_STAGES; i++)
			m_b_coeffs[i] = b[i];
		for (int i=0; i<T_N_INTERNAL_STAGES*T_N_INTERNAL_STAGES; i++)
			m_a_coeffs[i] = a[i];
		return true;
	}

	void
	setTruncatureOrder ( )
	{ m_order_q = std::max(0,T_N_INTERNAL_STAGES-2); }

	/**
	 * Computes the Runge-Kutta method.
	 */

	bool
	computeSolution (void)
	{
		bool success = false;

		T_LIE_ALGEBRA omega0 = T_LIE_ALGEBRA::Zero();
		T_LIE_ALGEBRA omega = omega0;
		T_LIE_ALGEBRA sol = omega0;
		T_LIE_ALGEBRA A;
		int i,j;

		for (i=0; i<T_N_INTERNAL_STAGES; i++) {

			// do something smarter if a(i,j) is zero
			omega = omega0;
			for (j=0; j<T_N_INTERNAL_STAGES; j++) {
				omega += this->m_k[j]*this->m_h*this->a_coeffs(i,j);
			}

			this->m_interface.computeA(A,omega.exp()*this->m_y0);
			// TODO: should you really overwrite m_k[i] ?
			this->m_k[i] = omega.computeDExpRInv(A,this->m_order_q);
			sol += this->m_b_coeffs[i]*this->m_h*this->m_k[i];
		}

		this->m_solution = sol;
		success = true;

		return success;
	}

	const T_Q
	reconstruct (const NOXVector<T_Q::DOF>& solution)
	{
		T_LIE_ALGEBRA w(solution);
		return T_Q(w.exp()*m_y0);
	}

	const T_Q
	reconstruct ()
	{ return this->reconstruct(this->m_solution.toVector()); }

	/*
	const T_Q
	reconstruct ()
	{
		T_LIE_ALGEBRA w = T_LIE_ALGEBRA::Zero();
		for (int i=0; i<T_N_INTERNAL_STAGES; i++)
			w += m_h*m_b_coeffs[i]*m_k[i];
		return T_Q(w.exp()*m_y0);
	}
	*/

	bool
	computeA (T_LIE_ALGEBRA& A, const T_Q& y)
	{ return m_interface.computeA(A,y); }

	bool
	computeJacobianA (std::vector<T_LIE_ALGEBRA>& JA, const T_Q& y)
	{ return m_interface.computeJacobianA(JA,y); }

};
} // namespace RKMK

#endif
