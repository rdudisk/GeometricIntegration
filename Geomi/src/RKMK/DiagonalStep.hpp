#ifndef DEF_RKMK_DIAGONALSTEP
#define DEF_RKMK_DIAGONALSTEP

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

namespace RKMK {

/**
 * This class inherits RKMK::Abstract::Step and implements a \f$p\f$ stages diagonal implicit RKMK step,
 * where \f$p\f$ is the template parameter T_N_INTERNAL_STAGES.
 * A step is diagonal implicit iff the coefficients of the Butcher tableau check \f$ a_{i,j}=0,\,\forall i\leq j\f$.
 *
 * The root \f$k_i\f$ of the equation \f$k_i-f\left(h\sum_{j=1}^ia_{i,j}k_j\right)=0\f$ is computed for every
 * ascending \f$i\f$ from \f$1\f$ to \f$p\f$ using the NOX library.
 * The configuration of the coefficients allows to solve for \f$p\f$ sets of \f$q\f$ equations
 * instead of one set of \f$p\times q\f$ equations, thus significantly decreasing the cost of the solver.
 */

template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES>
class DiagonalStep : public Abstract::Step<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>
{
private:
	Teuchos::RCP<Teuchos::ParameterList>			m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests;
	Teuchos::RCP<NOXGroup<T_Q,T_N_INTERNAL_STAGES>> m_grp;
	Teuchos::RCP<NOX::Solver::Generic>				m_solver;

public:
	DiagonalStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface)
	:	Abstract::Step<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>(interface)
	{
		// m_internals is already initalized by base constructor, overriding
		DiagonalStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>* internals = new DiagonalStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>(interface);
		this->m_internals = internals;

		// NOX
		this->m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
		Teuchos::ParameterList& solverParameters = *m_solverParametersPtr;
		solverParameters.set("Nonlinear Solver","Tensor Based");

		Teuchos::ParameterList& directionParameters = solverParameters.sublist("Direction");
		directionParameters.set("Method","Tensor");

		Teuchos::ParameterList& globalStrategyParameters = solverParameters.sublist("Line Search");
		globalStrategyParameters.set("Method","Curvilinear");

		Teuchos::ParameterList& lineSearchParameters = globalStrategyParameters.sublist(globalStrategyParameters.get("Method","Curvilinear"));
		lineSearchParameters.set("Lambda Selection","Halving");
		lineSearchParameters.set("Max Iters",20);

		Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,NOX::StatusTest::NormF::Unscaled));
		Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(50));
		this->m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));


		DiagonalStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>* tmp = static_cast<DiagonalStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>*>(this->m_internals);
		this->m_grp = Teuchos::rcp(new NOXGroup<T_Q,T_N_INTERNAL_STAGES>(*tmp));
		this->m_solver = NOX::Solver::buildSolver(this->m_grp,this->m_statusTests,this->m_solverParametersPtr);
	}

	~DiagonalStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> ()
	{ }

	const T_Q
	makeStep (void)
	{
		// TODO: faire des tests pour voir si tout se passe bien
		T_Q* Y1 = new T_Q();
		bool success = true;
		bool verbose = false;

		for (int i=0; i<T_N_INTERNAL_STAGES; i++) {
			try {
				(static_cast<DiagonalStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>*>(this->m_internals))->currentStep(i);
				NOX::StatusTest::StatusType status = m_solver->solve();
				const NOXGroup<T_Q,1>& solnGrp = dynamic_cast<const NOXGroup<T_Q,1>&>(m_solver->getSolutionGroup());
				const NOXVector<T_Q::DOF>& solnVec = dynamic_cast<const NOXVector<T_Q::DOF>&>(solnGrp.getX());
				(static_cast<DiagonalStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>*>(this->m_internals))->setSolutionK(solnVec);

				if (status != NOX::StatusTest::Converged)
					success = false;
			} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
		}

		success &= this->m_internals->computeSolution();
		*Y1 = this->m_internals->reconstruct();

		return *Y1;
	}
};
} // namespace RKMK


/*
template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES>
class DIRKMKStepInternals : public Abstract::RKMKStepInternals<T_M,T_Q,T_LIE_ALGEBRA>
{
private:
	int m_fevals;
	T_LIE_ALGEBRA m_initialGuess;
	T_LIE_ALGEBRA m_solution;

public:
	ExplicitRKMKStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface)
	: m_interface(interface)
	{
		m_fevals = 0;
		m_initialGuess = T_LIE_ALGEBRA::Zero();
		m_solution = T_LIE_ALGEBRA::Zero();
	}

	const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>
	getInitialGuess ()
	{
		return NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>::Zero();
	}

	bool
	computeF (NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>& f, const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>& x)
	{
		bool success = false;

		T_LIE_ALGEBRA omega0 = T_LIE_ALGEBRA::Zero();
		T_LIE_ALGEBRA omega = omega0;
		T_LIE_ALGEBRA A;
		NOXVector<T_Q::DOF> vec;
		int i,j;

		for (i=0; i<T_N_INTERNAL_STAGES; i++) {
			for (j=0; j<T_Q::DOF; j++) {
				m_k[i][j] = x[i*T_N_INTERNAL_STAGES+j];
			}
		}
		for (i=0; i<T_N_INTERNAL_STAGES; i++) {

			omega = omega0;
			for (j=0; j<T_N_INTERNAL_STAGES; j++) {
				omega += m_k[j]*m_h*this->a_coeffs(i,j);
			}

			m_interface.computeA(A,omega.exp()*m_y0);
			
			vec = (omega.computeDExpRInv(A,m_order_q)-m_k[i]).toVector();
			for(j=0; j<T_Q::DOF; j++) {
				f[i*T_N_INTERNAL_STAGES+j] = vec[j];
			}
		}

		m_fevals++;
		return true;
	}
};

template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES>
class DiagonalRKMKStep : public Abstract::RKMKStep<T_M,T_Q,T_LIE_ALGEBRA>
{
private:
	Teuchos::RCP<Teuchos::ParameterList> m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo> m_statusTests;
	Teuchos::RCP<NOXGroup<T_Q,T_N_INTERNAL_STAGES>> m_grp;
	Teuchos::RCP<NOX::Solver::Generic> m_solver;

public:
	DiagonalRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface,
																 std::vector<double> a,
																 std::vector<double> b)
	: Abstract::RKMKStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>(interface, a, b)
	{
		m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
		Teuchos::ParameterList& solverParameters = *m_solverParametersPtr;
		solverParameters.set("Nonlinear Solver","Tensor Based");

		Teuchos::ParameterList& directionParameters = solverParameters.sublist("Direction");
		directionParameters.set("Method","Tensor");

		Teuchos::ParameterList& globalStrategyParameters = solverParameters.sublist("Line Search");
		globalStrategyParameters.set("Method","Curvilinear");

		Teuchos::ParameterList& lineSearchParameters = globalStrategyParameters.sublist(globalStrategyParameters.get("Method","Curvilinear"));
		lineSearchParameters.set("Lambda Selection","Halving");
		lineSearchParameters.set("Max Iters",20);

		Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,NOX::StatusTest::NormF::Unscaled));
		Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(50));
		m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));

		Teuchos::RCP<NOXGroup<T_Q,T_N_INTERNAL_STAGES>> m_grp = Teuchos::rcp(new NOXGroup<T_Q,T_N_INTERNAL_STAGES>(m_internals));

		Teuchos::RCP<NOX::Solver::Generic> m_solver = NOX::Solver::buildSolver(m_grp,m_statusTests,m_solverParametersPtr);
	}

	~DiagonalRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> ()
	{ }

	const T_Q&
	makeStep (void)
	{
		// TODO: faire des tests pour voir si tout se passe bien
		T_Q* Y1 = new T_Q();
		bool success;
		bool verbose = true;

		try {
			NOX::StatusTest::StatusType status = m_solver->solve();
			const NOXGroup<T_Q,T_N_INTERNAL_STAGES>& solnGrp = dynamic_cast<const NOXGroup<T_Q,T_N_INTERNAL_STAGES>&>(m_solver->getSolutionGroup());
			const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>& solnVec = dynamic_cast<const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>&>(solnGrp.getX());
			*Y1 = m_internals.reconstruct(solnVec);

			if (status == NOX::StatusTest::Converged) {
				success = true;
			} else {
				success = false;
			}
		} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

		return *Y1;
	}
};
*/





// TEMPLATES

/*
template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES = 1>
class GenericRKMKInternals : public Abstract::Step<T_Q,T_N_INTERNAL_STAGES>
{
private:
	// Number of calls to computeF
	int m_fevals;
	// Initial guess
	T_LIE_ALGEBRA m_initialGuess;
	// Correct answer
	T_LIE_ALGEBRA m_solution;
	// Y0
	T_Q m_y0;
	// time step
	T_M m_h;

	// The interface
	Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& m_interface;

	bool m_isExplicit;
	bool m_isDIRK;

	// RK coeffs
	double			m_a_coeffs	[T_N_INTERNAL_STAGES * T_N_INTERNAL_STAGES];
	double			m_b_coeffs	[T_N_INTERNAL_STAGES];
	T_LIE_ALGEBRA	m_k			[T_N_INTERNAL_STAGES];
	unsigned int	order_q;


public:
	GenericRKMKInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface)
	: m_interface(interface)
	{
		m_fevals = 0;
		m_initialGuess = T_LIE_ALGEBRA::Zero();
		m_solution = T_LIE_ALGEBRA::Zero();
		//for (int i=0; i<T_N_INTERNAL_STAGES; i++)
		//	m_k[i] = T_LIE_ALGEBRA();
	}

	double
	a_coeffs (int i, int j) const
	{
		return m_a_coeffs[i*T_N_INTERNAL_STAGES+j];
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
		this->setIsExplicit();
		return true;
	}

	void
	setIsExplicit ( )
	{
		int i,j;
		bool tmp_isZero;

		m_isExplicit = true;
		m_isDIRK = true;
		for (i=0; i<T_N_INTERNAL_STAGES; i++) {
			m_isExplicit &= isZero<double>(m_a_coeffs[i*T_N_INTERNAL_STAGES+i]);
			for (j=i+1; j<T_N_INTERNAL_STAGES; j++) {
				tmp_isZero = isZero<double>(m_a_coeffs[i*T_N_INTERNAL_STAGES+j]);
				m_isExplicit &= tmp_isZero;
				m_isDIRK &= tmp_isZero;
			}
		}	
	}

	bool
	isExplicit ( ) const
	{
		return m_isExplicit;
	}

	void
	setTruncatureOrder ( )
	{
		m_order_q = std::max(0,T_N_INTERNAL_STAGES-2);
	}

	void
	setData (T_M h_var, T_Q y0_var)
	{
		m_h = h_var;
		m_y0 = y0_var;
	}

	const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>
	getInitialGuess ()
	{
		return NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>::Zero();
	}

	const T_Q
	reconstruct (const NOXVector<T_Q::DOF>& w)
	{
		// TODO
		return T_Q(T_LIE_ALGEBRA(w).exp()*m_y0);
	}

	const T_Q
	reconstruct ()
	{
		return this->reconstruct(this->m_solution.toVector());
	}

	bool
	computeDirectSolution (void)
	{
		if (!m_isExplicit)
			return false;

		bool success = false;

		T_LIE_ALGEBRA omega0 = T_LIE_ALGEBRA::Zero();
		T_LIE_ALGEBRA omega = omega0;
		T_LIE_ALGEBRA sol = omega0;
		T_LIE_ALGEBRA A;
		int i,j;

		for (i=0; i<T_N_INTERNAL_STAGES; i++) {

			omega = omega0;
			for (j=0; j<i; j++) {
				omega += m_k[j]*m_h*this->a_coeffs(i,j);
			}

			m_interface.computeA(A,omega.exp()*m_y0);
			m_k[i] = omega.computeDExpRInv(A,m_order_q);

			sol += m_b_coeffs[i]*m_k[i]*m_h;
		}
		m_solution = sol;
		success = true;
		m_fevals++;

		return success;
	}

	bool
	computeF (NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>& f, const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>& x)
	{
		bool success = false;

		T_LIE_ALGEBRA omega0 = T_LIE_ALGEBRA::Zero();
		T_LIE_ALGEBRA omega = omega0;
		T_LIE_ALGEBRA A;
		NOXVector<T_Q::DOF> vec;
		int i,j;

		for (i=0; i<T_N_INTERNAL_STAGES; i++) {
			for (j=0; j<T_Q::DOF; j++) {
				m_k[i][j] = x[i*T_N_INTERNAL_STAGES+j];
			}
		}
		for (i=0; i<T_N_INTERNAL_STAGES; i++) {

			omega = omega0;
			for (j=0; j<T_N_INTERNAL_STAGES; j++) {
				omega += m_k[j]*m_h*this->a_coeffs(i,j);
			}

			m_interface.computeA(A,omega.exp()*m_y0);
			
			vec = (omega.computeDExpRInv(A,m_order_q)-m_k[i]).toVector();
			for(j=0; j<T_Q::DOF; j++) {
				f[i*T_N_INTERNAL_STAGES+j] = vec[j];
			}
		}

		m_fevals++;
		return true;
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
	computeJacobianDIRK (Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>& J, const NOXVector<T_Q::DOF>& x, const int step_index)
	{
		m_k[step_index] = Algebra(x);

		int i,j,p;

		// Setting up variables
		
		double a = a_coeffs(step_index,step_index);
		// TODO : Check if a is zero and act accordingly
		
		T_LIE_ALGEBRA Gen;
		
		T_LIE_ALGEBRA w = T_LIE_ALGEBRA::Zero();
		for (j=0; j<=step_index; j++) {
			w += a_coeffs(step_index,i);
		}
		w = m_h*w;

		NOXVector<T_Q::DOF> Y = w.exp()*m_Y0;

		T_LIE_ALGEBRA A;
		computeA(A,Y);

		std::vector<T_LIE_ALGEBRA>& JA;
		computeJacobianA(JA,Y);
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
			dAdk = T_LIE_ALGEBRA(m_h*a*MatJA*w.partialExp(p)*m_Y0);
			dfdk = dAdk;
			// d(ad^k_w(A))/dk
			if (m_order_q > 0) {
				daddk = dAdk;
				ad = A;
				factorial = 1.0;
				for (i=0; i<m_order_q; i++) {
					factorial = factorial*(i+1);
					daddk = m_h*a*Gen.bracket(ad) + w.bracket(daddk);
					if (i%2==0) { // odd Bernoulli are 0, so no need for the update
						dfdk += (BERNOULLI_NUMBERS[i]/factorial)*daddk;
					}
					ad = w.bracket(ad);
				}
			}
			// col ? row ?
			J.col(p) = (dfdk-Gen).toVector();
		}
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF*T_N_INTERNAL_STAGES,T_Q::DOF*T_N_INTERNAL_STAGES>& J, const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>& x)
	{
		// TODO
		int i,j;
		return true;
	}
};
*/

/*
template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES = 1>
class GenericRKMKStep //: public Abstract::Step<T_Q>
{
private:
	// The interface
	Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& m_interface;

	GenericRKMKInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> m_internals;
	
	Teuchos::RCP<Teuchos::ParameterList> m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo> m_statusTests;
	Teuchos::RCP<NOXGroup<T_Q,T_N_INTERNAL_STAGES>> m_grp;
	Teuchos::RCP<NOX::Solver::Generic> m_solver;

public:
	GenericRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface)
	:	m_interface(interface),
		m_internals(GenericRKMKInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>(m_interface))
	{

		m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
		Teuchos::ParameterList& solverParameters = *m_solverParametersPtr;
		solverParameters.set("Nonlinear Solver","Tensor Based");

		Teuchos::ParameterList& directionParameters = solverParameters.sublist("Direction");
		directionParameters.set("Method","Tensor");

		Teuchos::ParameterList& globalStrategyParameters = solverParameters.sublist("Line Search");
		globalStrategyParameters.set("Method","Curvilinear");

		Teuchos::ParameterList& lineSearchParameters = globalStrategyParameters.sublist(globalStrategyParameters.get("Method","Curvilinear"));
		lineSearchParameters.set("Lambda Selection","Halving");
		lineSearchParameters.set("Max Iters",20);

		Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,NOX::StatusTest::NormF::Unscaled));
		Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(50));
		m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));

		Teuchos::RCP<NOXGroup<T_Q,T_N_INTERNAL_STAGES>> m_grp = Teuchos::rcp(new NOXGroup<T_Q,T_N_INTERNAL_STAGES>(m_internals));

		Teuchos::RCP<NOX::Solver::Generic> m_solver = NOX::Solver::buildSolver(m_grp,m_statusTests,m_solverParametersPtr);
	}

	~GenericRKMKStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> ()
	{ }

	void
	setIsExplicit ( )
	{
		m_internals.setIsExplicit();
	}

	bool
	isExplicit ( ) const
	{
		return m_internals.isExplicit();
	}

	void
	setData (T_M h_var, T_Q y0_var)
	{
		m_internals.setData(h_var,y0_var);
	}

	bool
	setCoeffs (std::vector<double> a, std::vector<double> b)
	{
		return m_internals.setCoeffs(a,b);
	}

	const T_Q&
	makeStep (void)
	{
		// TODO: faire des tests pour voir si tout se passe bien
		T_Q* Y1 = new T_Q();
		bool success;
		bool verbose = true;

		if (isExplicit()) {
			success = m_internals.computeDirectSolution();
			*Y1 = m_internals.reconstruct();
		}
		else {
			try {
				NOX::StatusTest::StatusType status = m_solver->solve();
				const NOXGroup<T_Q,T_N_INTERNAL_STAGES>& solnGrp = dynamic_cast<const NOXGroup<T_Q,T_N_INTERNAL_STAGES>&>(m_solver->getSolutionGroup());
				const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>& solnVec = dynamic_cast<const NOXVector<T_Q::DOF*T_N_INTERNAL_STAGES>&>(solnGrp.getX());
				*Y1 = m_internals.reconstruct(solnVec);

				if (status == NOX::StatusTest::Converged) {
					success = true;
				} else {
					success = false;
				}
			} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
		}

		return *Y1;
	}
};
*/

#endif
