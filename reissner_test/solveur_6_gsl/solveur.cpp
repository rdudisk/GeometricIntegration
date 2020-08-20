#include <iostream>
#include <cmath>

#include "classes.hpp"

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

#include <pugixml.hpp>

int
main (int argc, char* argv[])
{
	/* Setting up the problem ************************************************/
	

	pugi::xml_document doc;
	if (!doc.load_file("config.xml")) {
		std::cout << "Unable to open config.xml\nExiting" << std::endl;
		return -1;
	}

	double diameter, rho, young, poisson, area;

	pugi::xml_node xml_params = doc.child("config").child("params");
	young = std::stod(xml_params.child("young").child_value());
	poisson = std::stod(xml_params.child("poisson").child_value());
	rho = std::stod(xml_params.child("rho").child_value());
	diameter = std::stod(xml_params.child("diameter").child_value());
	area = diameter*diameter;

	double l, h;
	int n_space_steps, n_time_steps;

	pugi::xml_node xml_integration = doc.child("config").child("integration");
	l = std::stod(xml_integration.child("space-step").child_value());
	h = std::stod(xml_integration.child("time-step").child_value());
	n_space_steps = std::stoi(xml_integration.child("n-space-step").child_value());
	n_time_steps = std::stoi(xml_integration.child("n-time-step").child_value());

	RigidBody problem;

	problem.setInertia(area,rho);
	problem.setConstraint(area,young,poisson);
	problem.setSize(n_time_steps,n_space_steps);
	problem.baselinstep(0.0,h,0.0,l);

	/* Setup Stepper **********************************************************/

	/*
	Teuchos::RCP<Teuchos::ParameterList>			m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests;
	Teuchos::RCP<NOXGroup<NOXVector<3>,1>>			m_grp;
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
	m_grp = Teuchos::rcp(new NOXGroup<NOXVector<3>,1>(*solveme));
	m_solver = NOX::Solver::buildSolver(m_grp,m_statusTests,m_solverParametersPtr);
	*/

	/* Integration loop *******************************************************/

	int i, j;
	Group g0, g1;
	Algebra x0, x1, e0, e1;
	Vec6 mu0, mu1, sig0, sig1;
	Vec3 OM;
	bool verbose;
	int success;
	
	// Initialisation

	Eigen::Matrix<double,6,1> v_e0, v_e1, E4;
	v_e0 << 1.0, 1.0, 1.5, 1.0, 0.0, 0.0;
	v_e1 << 1.005, 1.004, 1.52, 1.0, 0.0, -0.01;
	E4 << 0, 0, 0, 1, 0, 0;
	e0 = Algebra(v_e0);
	e1 = Algebra(v_e1);

	g0 = Group::Identity();
	problem.pos(0,0,g0);
	g0.translation(2,h);
	problem.pos(1,0,g0);

	for (j=0; j<n_space_steps; j++) {
		if (j>0) {
			problem.pos(0,j,problem.pos(0,j-1)*Cay::eval(l*e0));
			problem.pos(1,j,problem.pos(1,j-1)*Cay::eval(l*e1));
		}
		x0 = (1.0/h)*Cay::inv(problem.pos(0,j).inverse()*problem.pos(1,j));
		mu0 = Cay::inv_right_diff_star(h*x0,Algebra(problem.Inertia()*x0.vector())).vector();
		//std::cout << "g_0^-1g_1" << std::endl << (problem.pos(0,j).inverse()*problem.pos(1,j)).matrix() << std::endl;
		//std::cout << "x0:" << std::endl << x0 << std::endl;
		//std::cout << "mu0:" << std::endl << mu0 << std::endl;
		problem.vel_time(0,j,x0);
		problem.mom_time(0,j,mu0);
		if (j<n_space_steps-1) {
			problem.vel_space(0,j,e0);
			problem.vel_space(1,j,e1);
		}
		else {
			problem.vel_space(0,j,Algebra(E4));
			problem.vel_space(1,j,Algebra(E4));
		}
	}

	// Integration
	
	int solve_err_nb = 0, max_err_nb = 20;
	Algebra* x_sol = new Algebra;
	
	for (i=1; i<n_time_steps-1; i++) {
		for (j=0; j<n_space_steps; j++) {
			if (j<n_space_steps-1) e1 = (1.0/l)*Cay::inv(problem.pos(i,j).inverse()*problem.pos(i,j+1));
			else e1 = Algebra(E4);
			problem.vel_space(i,j,e1);
			sig1 = (-1.0)*Cay::inv_right_diff_star(l*e1,Algebra(problem.Constraint()*(e1.vector()-E4))).vector();
			problem.mom_space(i,j,sig1);

			x0 = problem.vel_time(i-1,j);
			mu0 = problem.mom_time(i-1,j);

			mu1 = Algebra::static_Ad_star(Cay::eval(h*x0),Algebra(mu0)).vector() - (h/l)*sig1;

			if (j>0) {
				e0 = problem.vel_space(i,j-1);
				sig0 = problem.mom_space(i,j-1);
				mu1 += (h/l)*Algebra::static_Ad_star(Cay::eval(l*e0),Algebra(sig0)).vector();
			}

			problem.mom_time(i,j,mu1);

			std::cout << "Solving node " << i << ", " << j << "... ";
			success = solve_speed(mu1,h,x0,problem.Inertia(),x_sol);
			solve_err_nb += success;
			if (!success) std::cout << "success!" << std::endl;
			else std::cout << "fail!" << std::endl;

			x1 = *x_sol;
			problem.vel_time(i,j,x1);
			/*
			//if (i>1)
				//solveme->setOM((h/2.0)*(problem.vel_time(i-2,j).rot()),(h/2.0)*(problem.vel_time(i-1,j).rot()));
			//else
				solveme->setOM((h/2.0)*(problem.vel_time(i-1,j).rot()),(h/2.0)*(problem.vel_time(i-1,j).rot()));
			*/

			problem.vel_time(i,j,x1);
			problem.pos(i+1,j,problem.pos(i,j)*Cay::eval(h*x1));
		}
		// TODO: ceci n'est pas la bonne méthode, on fait ça en attendant mieux
		//e1 = Algebra(E4);
		//problem.pos(i+1,n_space_steps-1,problem.pos(i+1,n_space_steps-2)*(l*e1).cay());
		//x1 = (1.0/h)*Algebra::cay_inv(problem.pos(i,n_space_steps-1).inverse()*problem.pos(i+1,n_space_steps-1));
		//problem.vel_time(i+1,n_space_steps-1,x1);
		//mu1 = (h*x1).dCayRInv().transpose()*problem.Inertia().asDiagonal()*x1.toVector();

		/*
		if (success) {
			std::cout << "** ARRET **" << std::endl;
			break;
		}*/
	}

	/* Output results ********************************************************/

	problem.writeCSVFile("results.csv",true);

	//std::cout << "End result: TEST PASSED" << std::endl;

	return 0;
}

