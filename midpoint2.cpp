#include <iostream>
#include <cmath>

#include "include/RKMK"

#include "include/Space/SO3Group.hpp"
#include "include/Space/SO3Algebra.hpp"

#include <Eigen/Dense>
#include <Eigen/Geometry>

typedef Lie::SO3::Algebra<double> Algebra;
typedef Lie::SO3::Group<double> Group;

typedef NOXVector<3> Q;
typedef NOXVector<3> NOXV;
typedef NOXGroup<Q> NOXG;

class MyProblem : public RKMK::Abstract::Problem<double,Q,Algebra>
{
	private:
		Eigen::Matrix<double,3,1> m_Iinv;
	
	public:
		MyProblem()
		{
			m_Iinv << 1, 1, 1;
		}

		bool
		computeA (Algebra& A, const Q& x)
		{
			A = Algebra((-m_Iinv).asDiagonal()*x);
			return true;
		}

		bool
		computeJacobianA (std::vector<Algebra>& JA, const Q& x)
		{
			Eigen::Matrix<double,3,3> M = (-m_Iinv).asDiagonal();
			for (int i=0; i<3; i++)
				JA.push_back(Algebra(M.col(i)));
			return true;
		}
	
		Eigen::Matrix<double,3,1>&
		Iinv ()
		{
			return m_Iinv;
		}

		void
		Iinv (Eigen::Matrix<double,3,1> val)
		{
			m_Iinv = val;
		}
};

/*
typedef LieMidpointStep<double,Q,Algebra> LieMidStep;

class GenericIntegrator
{
private:
	Teuchos::RCP<Teuchos::ParameterList> m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo> m_statusTests;
	MyProblem& m_problem;
	
public:
	GenericIntegrator (MyProblem& problem)
	: m_problem(problem)
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
	}

	~GenericIntegrator ()
	{ }

	bool
	solve (void)
	{
		int i,j;
		bool success = false;
		bool verbose = true;
		double h;
		int n_steps = m_problem.size();
		Q Y0, Y1, tmp;

		for (i=0; i<n_steps-1; i++) {
			try {
				Y0 = m_problem.pos(i);
				h = m_problem.base(i+1)-m_problem.base(i);

				LieMidStep lieMidpointStep(m_problem,h,Y0);
				Teuchos::RCP<NOXG> grp = Teuchos::rcp(new NOXG(lieMidpointStep));

				Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,m_statusTests,m_solverParametersPtr);

				NOX::StatusTest::StatusType status = solver->solve();

				const NOXG& solnGrp = dynamic_cast<const NOXG&>(solver->getSolutionGroup());
				const NOXV& solnVec = dynamic_cast<const NOXV&>(solnGrp.getX());

				Y1 = lieMidpointStep.reconstruct(solnVec);
				m_problem.pos(i+1,Y1);

				if (status == NOX::StatusTest::Converged) {
					success = true;
				} else {
					success = false;
					break;
				}
			} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
		}
		return success;
	}

};
*/

////typedef GenericRKMKStep<double,Q,Algebra,4> RKMK4Step;
//typedef ExplicitRKMKStep<double,Q,Algebra,4> RKMK4Step;
//typedef ExplicitRKMKIntegrator<double,Q,Algebra> Integrator;

int
main (int argc, char* argv[])
{
	double h = 0.1;
	int n_steps = 4;
	MyProblem myProblem;
	myProblem.baselinstep(0.0,h,n_steps);

	Eigen::Matrix<double,3,1> Iinv(3.0/2.0,1.0,0.5);
	myProblem.Iinv(Iinv);

	Eigen::Matrix<double,3,1> pos(cos(M_PI/3.0), 0.0, sin(M_PI/3.0));
	myProblem.pos(0,pos);

	RKMK::Abstract::Integrator& integrator = RKMK::Factory<double,Q,Algebra>::createIntegrator(myProblem,"Implicit Euler");

	bool success = integrator.integrate();

	myProblem.write2csv("res_midpoint2.csv");

	return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

