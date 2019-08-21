#ifndef DEF_RKMK_DIAGONALSTEPINTERNALS
#define DEF_RKMK_DIAGONALSTEPINTERNALS

//#include "include/Step/AbstractStep.hpp"
//#include "RKMK_Abstract_Step.hpp"
//#include "../Common/Common_NOXVector.hpp"
//#include "../Common/Common_NOXGroup.hpp"

namespace RKMK {

template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES>
class DiagonalStepInternals : public StepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>, public ::Abstract::NOXStep<T_Q,1>
{
private:
	int m_fevals;
	T_LIE_ALGEBRA m_initialGuess;
	T_LIE_ALGEBRA m_solution;
	int m_currentStep;

public:
	DiagonalStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface)
	:	StepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>(interface)
	{
		m_fevals = 0;
		m_initialGuess = T_LIE_ALGEBRA::Zero();
		m_solution = T_LIE_ALGEBRA::Zero();
		m_currentStep = 0;
	}

	/*
	DIRKMKStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (const DIRKMKStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>& other)
	:	Abstract::RKMKStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>(interface)
	{
		m_fevals = other.m_fevals;
		m_initialGuess = other.m_initialGuess;
		m_solution = other.m_solution;
		m_currentStep = other.m_currentStep;
	}
	*/


	const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>
	getInitialGuess ()
	{ return NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>::Zero(); }

	int
	currentStep ()
	{ return m_currentStep; }

	// TODO : maybe some small checks ?
	void
	currentStep (int c)
	{ m_currentStep = c; }

	void
	setSolutionK (const NOXVector<T_Q::DOF>& x)
	{ this->m_k[m_currentStep] = T_LIE_ALGEBRA(x); }

	bool
	computeF (NOXVector<T_Q::DOF>& f, const NOXVector<T_Q::DOF>& x)
	{
		bool success = false;

		this->m_k[m_currentStep] = T_LIE_ALGEBRA(x);

		T_LIE_ALGEBRA omega0 = T_LIE_ALGEBRA::Zero();
		T_LIE_ALGEBRA omega = omega0;
		T_LIE_ALGEBRA res = omega0;
		T_LIE_ALGEBRA A;
		int i,j;

		omega = omega0;
		for (j=0; j<=m_currentStep; j++) {
			omega += this->m_k[j]*this->a_coeffs(m_currentStep,j);
		}
		omega = this->m_h*omega;
		this->m_interface.computeA(A,omega.exp()*this->m_y0);
		res = omega.computeDExpRInv(A,this->m_order_q) - this->m_k[m_currentStep];

		f = res.toVector();

		m_fevals++;
		return true;
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>& J, const NOXVector<T_Q::DOF>& x)
	{
		this->m_k[m_currentStep] = T_LIE_ALGEBRA(x);

		int i,j,p;
		double factorial;

		// Setting up variables
		
		double a = this->a_coeffs(m_currentStep,m_currentStep);
		// TODO : Check if a is zero and act accordingly
		
		T_LIE_ALGEBRA Gen;
		
		T_LIE_ALGEBRA w = T_LIE_ALGEBRA::Zero();
		for (j=0; j<=m_currentStep; j++) {
			w += this->a_coeffs(m_currentStep,j)*this->m_k[j];
		}
		w = this->m_h*w;

		NOXVector<T_Q::DOF> Y = w.exp()*this->m_y0;

		T_LIE_ALGEBRA A;
		this->computeA(A,Y);

		std::vector<T_LIE_ALGEBRA> JA;
		this->computeJacobianA(JA,Y);
		Eigen::Matrix<double,T_Q::DOF,T_Q::DOF> MatJA;
		for (i=0; i<T_Q::DOF; i++) {
			// TODO : row ? col ?
			MatJA.row(i) = (JA[i]).toVector();
		}

		T_LIE_ALGEBRA dAdk;
		T_LIE_ALGEBRA dfdk;
		T_LIE_ALGEBRA daddk;
		T_LIE_ALGEBRA ad;

		// Computing dF/dk_p

		for (p=0; p<T_Q::DOF; p++) {
			Gen = T_LIE_ALGEBRA::Generator(p);
			dAdk = T_LIE_ALGEBRA(this->m_h*a*MatJA*w.partialExp(p)*this->m_y0);
			dfdk = dAdk;
			// d(ad^k_w(A))/dk
			if (this->m_order_q > 0) {
				daddk = dAdk;
				ad = A;
				factorial = 1.0;
				for (i=0; i<this->m_order_q; i++) {
					factorial = factorial*(i+1);
					daddk = this->m_h*a*Gen.bracket(ad) + w.bracket(daddk);
					if (i%2==0) { // odd Bernoulli are 0, so no need for the update
						dfdk += (BERNOULLI_NUMBERS[i]/factorial)*daddk;
					}
					ad = w.bracket(ad);
				}
			}
			// col ? row ?
			J.col(p) = (dfdk-Gen).toVector();
		}

		return true;
	}
};
}

#endif
