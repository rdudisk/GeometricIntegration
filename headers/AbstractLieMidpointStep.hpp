#ifndef DEF_ABSTRACT_LIE_MIDPOINT_SOLVER
#define DEF_ABSTRACT_LIE_MIDPOINT_SOLVER

class AbstractLieMidpointSolver : public NOX::LAPACK::Interface
{
private:
	// Number of calls to computeF
	int fevals;
	// Initial guess
	NOX::LAPACK::Vector initialGuess;
	// Correct answer
	NOX::LAPACK::Vector solution;
	// Y0
	Algebra y0;
	// Iinv
	Eigen::Matrix<double,3,1> Iinv;
	// time step
	double h;

public:

	LieMidpointSolver (double h_var, Algebra y0_var, Eigen::Matrix<double,3,1> Iinv_var) :
		h(h_var),
		y0(y0_var),
		Iinv(Iinv_var),
		initialGuess(3),
		solution(3)
	{
		fevals = 0;

		for (int i=0; i<3; i++) {
			initialGuess(i) = 0;
			solution(i) = 1;
		}

		std::cout << "Rigid body problem: y0 = " << y0 << std::endl;
	}

	~LieMidpointSolver ()
	{
		std::cout << "Function evaluations: " << fevals << std::endl;
	}

	const NOX::LAPACK::Vector&
	getInitialGuess (void)
	{
		return initialGuess;
	}

	const NOX::LAPACK::Vector&
	getSolution (void)
	{
		return solution;
	}

	virtual bool
	computeA (AbstractAlgebra& A, const NOX::Abstract::Vector& x) = 0;

	bool
	computeF (NOX::Abstract::Vector& f, const NOX::Abstract::Vector& x)
	{
		int i;
		AbstractAlgebra w0(x);
		AbstractAlgebra w1(w*0.5);
		AbstractAlgebra A;
		computeA(A,w1.exp()*y0);
		AbstractAlgebra res(h*A-w);

		for (i=0; i<3; i++) {
			f(i) = (res.toVector())(i);
		}
		
		fevals++;
		return true;
	}

	virtual bool
	computeJacobianA (std::vector<AbstractAlgebra> JA, const NOX::Abstract::Vector& x) = 0;

	/*
	 * Voir calculs cahier III p.80
	 */
	bool
	computeJacobian (NOX::LAPACK::Matrix<double>& J, const NOX::LAPACK::Vector& x)
	{
		int i,j;
		AbstractAlgebra w(x);
		AbstractAlgebra w1(w*0.5);

		double nm = chi.norm();
		Eigen::Matrix<double,3,3> K = (chi.normalized()).toRotationMatrix();
		Eigen::Matrix<double,3,3> dexp;
		Algebra dA, df;
		Eigen::Matrix<double,3,1> dF;
		double c = cos(nm), s = sin(nm);

		for (j=0;j<3;j++) {
			Eigen::Matrix<double,3,3> M = Algebra::GeneratorMatrix(j);
			dexp = (isZero(nm)) ? M : (c-(s/nm))*(chi[j]/nm)*K + (s/nm)*M + (s+(1.0-c)*(2.0/nm))*(chi[j]/nm)*K*K + ((1-c)/nm)*(M*K+K*M);
			dA = Algebra(-0.5*II*dexp*(y0.toVector()));
			df = h*dA-Algebra::Generator(j);
			dF = df.toVector();
			for (i=0;i<3;i++) {
				J(i,j) = dF(i);
			}
		}

		return true;
	}

};

#endif
