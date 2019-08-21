#ifndef DEF_LIE_MIDPOINT_STEP
#define DEF_LIE_MIDPOINT_STEP

/*
#include "include/Problem/ProblemInterface.hpp"
#include "include/Common/NOXVector.hpp"
#include "include/Step/AbstractStep.hpp"

bool isZero (double);

template <typename T_M, typename T_Q, typename T_LIE_ALGEBRA>
class LieMidpointStep : public Abstract::Step<T_Q>
{
private:
	// Number of calls to computeF
	int fevals;
	// Initial guess
	T_LIE_ALGEBRA initialGuess;
	// Correct answer
	T_LIE_ALGEBRA solution;
	// Y0
	T_Q y0;
	// time step
	T_M h;
	// The interface
	Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& m_interface;

public:

	LieMidpointStep<T_M,T_Q,T_LIE_ALGEBRA>
		(	Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface,
			T_M h_var,
			T_Q y0_var) :
		m_interface(interface),
		h(h_var),
		y0(y0_var)
	{
		fevals = 0;

		initialGuess = T_LIE_ALGEBRA::Zero();
		solution = T_LIE_ALGEBRA::Zero();
	}

	~LieMidpointStep<T_M,T_Q,T_LIE_ALGEBRA> ()
	{ }

	void
	setData (T_M h_var, T_Q y0_var)
	{
		h = h_var;
		y0 = y0_var;
	}

	T_Q
	reconstruct (const NOXVector<T_Q::DOF>& w)
	{
		return T_Q(T_LIE_ALGEBRA(w).exp()*y0);
	}

	const NOXVector<T_Q::DOF>
	getInitialGuess (void)
	{
		return initialGuess.toNOXVector();
	}

	const T_LIE_ALGEBRA&
	getSolution (void)
	{
		return solution;
	}

	bool
	computeA (T_LIE_ALGEBRA& A, const T_Q& y)
	{
		return m_interface.computeA(A,y);
	}

	bool
	computeJacobianA (std::vector<T_LIE_ALGEBRA>& JA, const T_Q& y)
	{
		return m_interface.computeJacobianA(JA,y);
	}

	bool
	computeF (NOXVector<T_Q::DOF>& f, const NOXVector<T_Q::DOF>& x)
	{
		int i;
		T_LIE_ALGEBRA w0(x);
		T_LIE_ALGEBRA w1(w0*0.5);
		T_LIE_ALGEBRA A;
		m_interface.computeA(A,w1.exp()*y0);
		T_LIE_ALGEBRA res(h*A-w1);

		f = res.toVector();
		
		fevals++;
		return true;
	}

	// Voir calculs cahier III p.80 et surtout p.95
	bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>& J, const NOXVector<T_Q::DOF>& x)
	{
		int i,j;
		T_LIE_ALGEBRA w(x);
		T_LIE_ALGEBRA w1(w*0.5);

		T_Q Y = w1.exp()*y0;

		std::vector<T_LIE_ALGEBRA> JA;
		this->computeJacobianA(JA,Y);

		Eigen::Matrix<double,3,3> partialExp;
		Eigen::Matrix<double,3,1> partialF;

		for (i=0; i<T_Q::DOF; i++) {
			partialF = -T_LIE_ALGEBRA::Generator(i).toVector();
			partialExp = w1.partialExp(i);
			for (j=0; j<T_Q::DOF; j++) {
				partialF += 2.0*h*((JA[j]).toVector()*(partialExp.row(j)*y0));
			}
			J.col(i) = partialF;
		}

		return true;
	}

};

*/
#endif
