#include <iostream>
#include <cmath>

#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_LAPACK_Group.H"
#include "Teuchos_StandardCatchMacros.hpp"

#include "headers/SO3Group.hpp"
#include "headers/SO3Algebra.hpp"
#include "headers/DiscSystem.hpp"

#include <Eigen/Dense>
#include <Eigen/Geometry>

// Voir https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/float_comparison.html
bool
isZero (double d)
{
	return (abs(d)<std::numeric_limits<double>::epsilon());
}

typedef Lie::SO3::Algebra<double> Algebra;
typedef Lie::SO3::Group<double> Group;

Algebra
A (const Eigen::Matrix<double,3,1>& y, const Eigen::Matrix<double,3,1> Iinv)
{
	return Algebra((-Iinv).asDiagonal()*y);
}

class LieMidpointSolver : public NOX::LAPACK::Interface
{
private:
	// Numer of calls to computeF
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

	bool
	computeF (NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector& x)
	{
		int i;
		Eigen::Matrix<double,3,1> v;
		
		for (i=0; i<3; i++) {
			v(i) = x(i);
		}

		Algebra w(v);
		Algebra w2(w*0.5);
		Algebra res;
		Eigen::Matrix<double,3,3> II(Iinv.asDiagonal());
		Algebra A(-II*((w2.exp()).transformVector(y0.toVector())));

		res = h*A-w;
		for (i=0; i<3; i++) {
			f(i) = (res.toVector())(i);
		}
		
		fevals++;
		return true;
	}

	/*
	 * Voir calculs cahier III pp.74-77
	 */
	bool
	computeJacobian (NOX::LAPACK::Matrix<double>& J, const NOX::LAPACK::Vector& x)
	{
		int i,j;
		Eigen::Matrix<double,3,1> v;
		
		for (i=0; i<3; i++) {
			v(i) = x(i);
		}

		Algebra w(v);
		Algebra chi(v*0.5);
		Eigen::Matrix<double,3,3> II(Iinv.asDiagonal());

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

class MidpointIntegrator
{
private:
	Teuchos::RCP<Teuchos::ParameterList> m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo> m_statusTests;
	DiscSyst<double,Algebra>& m_system;
	
public:
	MidpointIntegrator (DiscSyst<double,Algebra>& system)
	: m_system(system)
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

		/*Teuchos::ParameterList& printParams = solverParameters.sublist("Printing");
		  printParams.set("Output Precision",3);
		  printParams.set("Output Processor",0);
		  printParams.set("Output Information",
		  NOX::Utils::OuterIteration +
		  NOX::Utils::OuterIterationStatusTest +
		  NOX::Utils::InnerIteration +
		  NOX::Utils::Parameters + 
		  NOX::Utils::Details +
		  NOX::Utils::Warning);
		  NOX::Utils utils(printParams);
		  */

		Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-12,NOX::StatusTest::NormF::Unscaled));
		Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(50));
		m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));
	}

	~MidpointIntegrator ()
	{ }

	bool
	solve (void)
	{
		int i,j;
		bool success = true;
		double h;
		int n_steps = m_system.size();
		Algebra Y0, Y1, tmp;

		Eigen::Matrix<double,3,1> Iinv;
		Iinv << 3.0/2.0, 1.0, 1.0/2.0;

		for (i=0; i<n_steps-1; i++) {
			//try {
			Y0 = m_system.pos(i);
			h = m_system.base(i+1)-m_system.base(i);

			LieMidpointSolver lieMidpointSolver(h,Y0,Iinv);
			Teuchos::RCP<NOX::LAPACK::Group> grp = Teuchos::rcp(new NOX::LAPACK::Group(lieMidpointSolver));

			Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,m_statusTests,m_solverParametersPtr);

			NOX::StatusTest::StatusType status = solver->solve();

			const NOX::LAPACK::Group& solnGrp = dynamic_cast<const NOX::LAPACK::Group&>(solver->getSolutionGroup());
			const NOX::LAPACK::Vector& solnVec = dynamic_cast<const NOX::LAPACK::Vector&>(solnGrp.getX());

			for (j=0; j<3; j++) {
				tmp[j] = solnVec(j);
			}

			Y1 = Algebra(Algebra(tmp).exp().transformVector(Y0.toVector()));
			m_system.pos(i+1,Y1);

			if (status == NOX::StatusTest::Converged) {
				success = true;
			} else {
				success = false;
				break;
			}
			//} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
		}
		return success;
	}

};

int
main (int argc, char* argv[])
{
	double h = 0.1;
	int n_steps = 4;
	DiscSyst<double,Algebra> ds;
	ds.baselinstep(0.0,h,n_steps);

	Eigen::Matrix<double,3,1> v;
	v << cos(M_PI/3.0), 0.0, sin(M_PI/3.0);
	//djs.pos(0,Algebra(v));
	ds.pos(0,Algebra(cos(M_PI/3.0), 0.0, sin(M_PI/3.0)));

	MidpointIntegrator integrator(ds);

	bool success = integrator.solve();

	ds.write2csv("res_midpoint.csv");

	return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
