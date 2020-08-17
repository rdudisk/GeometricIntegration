#include <iostream>
#include <cmath>

#include "classes.hpp"

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

int
main (int argc, char* argv[])
{
	/* Setting up the problem ************************************************/

	RigidBody problem;

	double diameter = 0.01;
	double area = diameter*diameter;
	double rho = 1.0e3;
	double young = 5.0e3;
	double poisson = 0.35;
	problem.setInertia(area,rho);
	problem.setConstraint(area,young,poisson);

	// linspace
	double l = 0.1;
	int n_space_steps = 10;
	double h = 0.0005; //l*problem.coeffCFL(young,poisson,rho,10.0);
	int n_time_steps = 300;
	problem.setSize(n_time_steps,n_space_steps);
	problem.baselinstep(0.0,h,0.0,l);

	/* Setup Stepper **********************************************************/

	Teuchos::RCP<Teuchos::ParameterList>			m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests;
	Teuchos::RCP<NOXGroup<NOXVector<6>,1>>			m_grp;
	Teuchos::RCP<NOX::Solver::Generic>				m_solver;

	SolveMe* solveme = new SolveMe(problem);

	m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
	Teuchos::ParameterList& solverParameters = *m_solverParametersPtr;

	solverParameters.set("Nonlinear Solver","Line Search Based");
	Teuchos::ParameterList& lineSearchParameters = solverParameters.sublist("Line Search");
	lineSearchParameters.set("Method","Full Step");


	if (true) { // true = silent
		Teuchos::ParameterList& printParams = solverParameters.sublist("Printing");
		printParams.set("Output Information",
			NOX::Utils::TestDetails +
			NOX::Utils::Error);
	}

	Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,NOX::StatusTest::NormF::Unscaled));
	Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(500));
	m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));

	// MidpointStepInternals<T_M,T_Q,T_TQ>* tmp = static_cast<MidpointStepInternals<T_M,T_Q,T_TQ>*>(this->m_internals);
	// m_grp = Teuchos::rcp(new NOXGroup<T_Q,1>(*tmp));
	m_grp = Teuchos::rcp(new NOXGroup<NOXVector<6>,1>(*solveme));
	m_solver = NOX::Solver::buildSolver(m_grp,m_statusTests,m_solverParametersPtr);

	/* Integration loop *******************************************************/

	int i, j;
	Group g0, g1;
	Algebra x0, x1, e0, e1;
	Vec6 mu0, mu1, sig0, sig1;
	Vec3 OM;
	bool success, verbose;
	
	// Initialisation

	Eigen::Matrix<double,6,1> v_e0, v_e1, E4;
	v_e0 << 1.0, 1.5, 1.0, 1.0, 0.0, 0.0;
	v_e1 << 1.004, 1.52, 1.005, 1.0, 0.0, -0.01;
	E4 << 0, 0, 0, 1, 0, 0;
	e0 = Algebra(v_e0);
	e1 = Algebra(v_e1);

	g0 = Group::Identity();
	problem.pos(0,0,g0);
	g0.trans(2,h);
	problem.pos(1,0,g0);

	for (j=0; j<n_space_steps; j++) {
		if (j>0) {
			problem.pos(0,j,problem.pos(0,j-1)*(l*e0).cay());
			problem.pos(1,j,problem.pos(1,j-1)*(l*e1).cay());
		}
		x0 = (1.0/h)*Algebra::cay_inv(problem.pos(0,j).inverse()*problem.pos(1,j));
		problem.vel_time(0,j,x0);
		problem.mom_time(0,j,(h*x0).dCayRInv().transpose()*problem.Inertia()*x0.toVector());
		if (j<n_space_steps-1) {
			problem.vel_space(0,j,e0);
			problem.vel_space(1,j,e1);
		}
		else {
			problem.vel_space(0,j,Algebra(E4));
			problem.vel_space(1,j,Algebra(E4));
		}
		std::cout << "*******************"
			<< "\nnode 0, " << j
			<< "\n\tpos " << problem.pos(0,j).matrix()
			<< "\n\tvel_time " << problem.vel_time(0,j)
			<< "\n\tvel_space " << problem.vel_space(0,j) << std::endl;
	}

	// Integration
	
	for (i=1; i<n_time_steps-1; i++) {
		for (j=0; j<n_space_steps; j++) {
			if (j<n_space_steps-1) e1 = (1.0/l)*Algebra::cay_inv(problem.pos(i,j).inverse()*problem.pos(i,j+1));
			else e1 = Algebra(E4);
			problem.vel_space(i,j,e1);
			sig1 = -((l*e1).dCayRInv().transpose()*problem.Constraint()*(e1.toVector()-E4));
			problem.mom_space(i,j,sig1);

			x0 = problem.vel_time(i-1,j);
			mu0 = problem.mom_time(i-1,j);

			mu1 = Algebra::static_Ad_star((h*x0).cay(),Algebra(mu0)).toVector() - (h/l)*sig1;

			if (j>0) {
				e0 = problem.vel_space(i,j-1);
				sig0 = problem.mom_space(i,j-1);
				mu1 += (h/l)*Algebra::static_Ad_star((l*e0).cay(),Algebra(sig0)).toVector();
			}

			problem.mom_time(i,j,mu1);

			solveme->setData(h,mu1);
			//if (i>1)
				//solveme->setOM((h/2.0)*(problem.vel_time(i-2,j).rot()),(h/2.0)*(problem.vel_time(i-1,j).rot()));
			//else
				//solveme->setOM((h/2.0)*(problem.vel_time(i-1,j).rot()),(h/2.0)*(problem.vel_time(i-1,j).rot()));
			solveme->setx0(x0.toVector());

			success = true;
			verbose = false;
			std::cout << "Solving node " << i << ", " << j << "... ";
			try {
				m_solver->reset(solveme->getInitialGuess());
				NOX::StatusTest::StatusType status = m_solver->solve();
				const NOXGroup<NOXVector<6>,1>& solnGrp = dynamic_cast<const NOXGroup<NOXVector<6>,1>&>(m_solver->getSolutionGroup());
				const NOXVector<6>& c_x = dynamic_cast<const NOXVector<6>&>(solnGrp.getX());

				if (status != NOX::StatusTest::Converged)
					success = false;

				x1 = c_x;
			} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
			if (success) std::cout << "success!" << std::endl;
			else std::cout << "fail!" << std::endl;

			std::cout << "*******************"
				<< "\nnode 0, " << j
				<< "\n\tpos " << problem.pos(0,j).matrix()
				<< "\n\tvel_time " << problem.vel_time(0,j)
				<< "\n\tvel_space " << problem.vel_space(0,j) << std::endl;
			// Update
			problem.vel_time(i,j,x1);
			problem.pos(i+1,j,problem.pos(i,j)*(h*x1).cay());
		}
		// TODO: ceci n'est pas la bonne méthode, on fait ça en attendant mieux
		//e1 = Algebra(E4);
		//problem.pos(i+1,n_space_steps-1,problem.pos(i+1,n_space_steps-2)*(l*e1).cay());
		//x1 = (1.0/h)*Algebra::cay_inv(problem.pos(i,n_space_steps-1).inverse()*problem.pos(i+1,n_space_steps-1));
		//problem.vel_time(i+1,n_space_steps-1,x1);
		//mu1 = (h*x1).dCayRInv().transpose()*problem.Inertia().asDiagonal()*x1.toVector();

		if (!success) {
			std::cout << "** ARRET **" << std::endl;
			break;
		}
	}

	/* Output results ********************************************************/

	problem.writeCSVFile("results.csv",true);

	//std::cout << "End result: TEST PASSED" << std::endl;

	return 0;
}

