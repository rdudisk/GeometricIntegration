#include <iostream>
#include <cmath>
#include <string>

#include "classes.hpp"

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

#include "pugixml.hpp"

int
main (int argc, char* argv[])
{
	/* Setting up the problem ************************************************/

	pugi::xml_document doc;
	if (!doc.load_file("config.xml")) {
		std::cout << "Unable to open config.xml\nExiting" << std::endl;
		return -1;
	}

	double radius, rho, young, poisson, k_factor, f_factor, pincement, l, h;
	int    n_space_steps, n_time_steps;

	pugi::xml_node xml_params = doc.child("config").child("params");
	young     = std::stod(xml_params.child("young").child_value());
	poisson   = std::stod(xml_params.child("poisson").child_value());
	rho       = std::stod(xml_params.child("rho").child_value());
	radius    = std::stod(xml_params.child("radius").child_value());
	k_factor  = std::stod(xml_params.child("kfactor").child_value());
	f_factor  = std::stod(xml_params.child("ffactor").child_value());
	pincement = std::stod(xml_params.child("pincement").child_value());

	pugi::xml_node xml_integration = doc.child("config").child("integration");
	l = std::stod(xml_integration.child("space-step").child_value());
	h = std::stod(xml_integration.child("time-step").child_value());
	n_space_steps = std::stoi(xml_integration.child("n-space-step").child_value());
	n_time_steps  = std::stoi(xml_integration.child("n-time-step").child_value());

	/* Checking up CFL condition requirement **********************************/

	double alpha = 0.1;
	double h_cfl=l*RigidBody::coeffCFL(young, poisson, rho, alpha);
	
	if (h>=h_cfl) {
		std::string response;
		std::cout << "Warning: CFL condition not verified (factor=" << alpha << "):" << std::endl
			      << "\th = " << h << std::endl
				  << "\tlimit = " << h_cfl << std::endl
				  << "Do you want to continue ? [yn] ";
		std::cin >> response;
		if (response!="y")
			return 0;
	}

	/* Instantiate Reissner beam problem **************************************/

	RigidBody problem;

	/*
	problem.setInertia(area,rho);
	problem.setConstraint(area,young,poisson);
	problem.setSize(n_time_steps,n_space_steps);
	problem.baselinstep(0.0,h,0.0,l);
	*/

	problem.setInertia(radius,rho);
	problem.setConstraint(radius,young,poisson);
	problem.setSize(n_space_steps);
	problem.setStepSize(h,l);
	problem.setCSV("results.csv");

	Displacement displacement(1.0/h,"displacement.csv");

	// TODO: clean that up
	double area = RigidBody::compute_area(radius);
	double F_constant = f_factor/(rho*area);
	double F_duration = 0.01;
	Eigen::Matrix<double,6,1> Fvector;
	Fvector << 0.0, 0.0, 0.0, 0.0, F_constant, 0.0;
	int indice_pincement = (int) (pincement*((float)n_space_steps));
	int indice_ecoute = indice_pincement;

	/* Setting up solver ******************************************************/

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

	Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-12,NOX::StatusTest::NormF::Scaled));
	Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(1000));
	m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));

	m_grp = Teuchos::rcp(new NOXGroup<NOXVector<6>,1>(*solveme));
	m_solver = NOX::Solver::buildSolver(m_grp,m_statusTests,m_solverParametersPtr);

	/* Solver for boudary j=M *************************************************/

	Teuchos::RCP<Teuchos::ParameterList>			m_solver2ParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests2;
	Teuchos::RCP<NOXGroup<NOXVector<6>,1>>			m_grp2;
	Teuchos::RCP<NOX::Solver::Generic>				m_solver2;

	SolveBoundary* solvebound = new SolveBoundary(problem);

	m_solver2ParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
	Teuchos::ParameterList& solver2Parameters = *m_solver2ParametersPtr;

	solver2Parameters.set("Nonlinear Solver","Line Search Based");
	Teuchos::ParameterList& lineSearchParameters2 = solver2Parameters.sublist("Line Search");
	lineSearchParameters2.set("Method","Full Step");

	if (true) { // true = silent
		Teuchos::ParameterList& printParams2 = solver2Parameters.sublist("Printing");
		printParams2.set("Output Information",
			NOX::Utils::TestDetails +
			NOX::Utils::Error);
	}

	Teuchos::RCP<NOX::StatusTest::NormF> statusTest2A = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-12,NOX::StatusTest::NormF::Scaled));
	Teuchos::RCP<NOX::StatusTest::MaxIters> statusTest2B = Teuchos::rcp(new NOX::StatusTest::MaxIters(1000));
	m_statusTests2 = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTest2A,statusTest2B));

	m_grp2 = Teuchos::rcp(new NOXGroup<NOXVector<6>,1>(*solvebound));
	m_solver2 = NOX::Solver::buildSolver(m_grp2,m_statusTests2,m_solver2ParametersPtr);

	/* Initialization *********************************************************/

	int i,j;
	bool verbose;
	bool success;
	double w_resample = (1.0/44100.0);
	
	Group   g0;
	Algebra x0,  x1,  e0,   e1;
	Vec6    mu0, mu1, sig0, sig1;
	Eigen::Matrix<double,6,1> E4;

	// Setting up a equilibrium reference beam with null speed
	
	g0   = Group::Identity();
	E4   << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
	e0   = Algebra(E4);
	x0   = Algebra::Zero();
	x1   = Algebra::Zero();
	x1[4]= 0.01*h;
	mu0  = Vec6::Zero();
	sig0 = Vec6::Zero();

	for (j=0; j<n_space_steps; j++) {
		g0.translation(0,j*l);
		problem.pos(0,j,g0);
		if (false) { //(j==3) {
			problem.vel_time(0,j,x1);
			problem.mom_time(0,j,Cay::inv_right_diff_star(h*x1,Algebra(problem.Inertia()*(x1.vector()))).vector());
		}
		else {
			problem.vel_time(0,j,x0);
			problem.mom_time(0,j,mu0);
		}
	}

	for (j=0; j<n_space_steps-1; j++) {
		e0   = (1.0/l)*Cay::inv(problem.pos(0,j).inverse()*problem.pos(0,j+1));
		sig0 = (-1.0)*Cay::inv_right_diff_star(l*e0,Algebra(problem.Constraint()*(e0.vector()-E4))).vector();
		problem.vel_space(0,j,e0);
		problem.mom_space(0,j,sig0);
	}
	
	displacement.insert(problem.pos(0,indice_ecoute).translationVector()[1]);


	// Updating position for i=1
	
	for (j=0; j<n_space_steps; j++) {
		problem.pos(1,j,problem.pos(0,j)*Cay::eval(h*problem.vel_time(0,j)));
	}
	
	problem.updateCSV(0,w_resample);

	/* Integration loop *******************************************************/

	// TODO: clean that up
	
	int solve_err_nb = 0, max_err_nb = 20;
	Algebra* x_sol = new Algebra;
	Eigen::Matrix<double,6,6> Kmatrix;
	Eigen::Matrix<double,6,1> Kvector;
	Kvector << 0.0, 0.0, 0.0, 1.0, 1.0, 1.0;
	Kmatrix = k_factor*Kvector.asDiagonal();
	double i_double;
	
	for (i=1; i<n_time_steps-1; i++) {
		for (j=0; j<n_space_steps; j++) {
			if (j<n_space_steps-1) {
				// Updating space velocity and momentum
				e1   = (1.0/l)*(Cay::inv(problem.pos(i,j).inverse()*problem.pos(i,j+1)));
				sig1 = (-1.0)*(Cay::inv_right_diff_star(l*e1,Algebra(problem.Constraint()*(e1.vector()-E4))).vector());
				problem.vel_space(i,j,e1);
				problem.mom_space(i,j,sig1);
			}
			//else
				//e1 = Algebra(E4);

			if (j==0) {
				// Updating time velocity and momentum
				x1  = Algebra::Zero();
				mu1 = Vec6::Zero();
				problem.vel_time(i,j,x1);
				problem.mom_time(i,j,mu1);

				// Updating position
				problem.pos(i+1,j,problem.pos(i,j));
			}
			else if (j==(n_space_steps-1)) {
				// solve for epsilon(i+1,M-1)
				solvebound->setData(l,
						problem.pos(i+1,j-1).rotationMatrix(),
						problem.pos(i+1,j-1).translationVector()-problem.pos(i,j).translationVector(),
						E4,
						problem.vel_space(i,j-1).vector());
				try {
					m_solver2->reset(solvebound->getInitialGuess());
					NOX::StatusTest::StatusType status = m_solver2->solve();
					const NOXGroup<NOXVector<6>,1>& solnGrp = dynamic_cast<const NOXGroup<NOXVector<6>,1>&>(m_solver2->getSolutionGroup());
					const NOXVector<6>& c_e = dynamic_cast<const NOXVector<6>&>(solnGrp.getX());

					if (status == NOX::StatusTest::Failed || status == NOX::StatusTest::Unconverged)
						success = false;
					else success = true;

					e1 = Algebra(c_e);
				} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
				if (!success) {
					std::cout << "Failed solving epsilon" << std::endl;
					return 0;
				}

				// Updating space velocity and momentum
				sig1 = (-1.0)*(Cay::inv_right_diff_star(l*e1,Algebra(problem.Constraint()*(e1.vector()-E4))).vector());
				problem.vel_space(i+1,j-1,e1);
				problem.mom_space(i+1,j-1,sig1);

				// Updating position
				problem.pos(i+1,j,problem.pos(i+1,j-1)*Cay::eval(l*e1));

				// Updating time velocity and momentum
				x1  = (1.0/h)*(Cay::inv(problem.pos(i,j).inverse()*problem.pos(i+1,j)));
				mu1 = Cay::inv_right_diff_star(h*x1,Algebra(problem.Inertia()*(x1.vector()))).vector();
				problem.vel_time(i,j,x1);
				problem.mom_time(i,j,mu1);
			}
			else {
				x0   = problem.vel_time(i-1,j);
				mu0  = problem.mom_time(i-1,j);
				e0   = problem.vel_space(i,j-1);
				sig0 = problem.mom_space(i,j-1);

				// Writing DEL equation
				mu1 = Algebra::static_Ad_star(Cay::eval(h*x0),Algebra(mu0)).vector();
				mu1 += (h/l)*(Algebra::static_Ad_star(Cay::eval(l*e0),Algebra(sig0)).vector() -sig1);
				// Adding friction forces
				mu1 -= h*l*Kmatrix*x0.vector();
				// Adding external control forces
				i_double = (double) i;
				if (j==indice_pincement && i_double*h < F_duration)
					mu1 += h*l*(i_double*h/F_duration)*Fvector;

				// Solving for x1
				solveme->setData(h,mu1,x0.vector());
				try {
					m_solver->reset(solveme->getInitialGuess());
					NOX::StatusTest::StatusType status = m_solver->solve();
					const NOXGroup<NOXVector<6>,1>& solnGrp = dynamic_cast<const NOXGroup<NOXVector<6>,1>&>(m_solver->getSolutionGroup());
					const NOXVector<6>& c_x = dynamic_cast<const NOXVector<6>&>(solnGrp.getX());

					if (status == NOX::StatusTest::Failed || status == NOX::StatusTest::Unconverged)
						success = false;
					else success = true;

					x1 = Algebra((2.0/h)*c_x);
				} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
				//if (success) {
					//std::cout << "Success with" << std::endl;
					//std::cout << "M:  " << 0.5*h*mu1 << std::endl;
					//std::cout << "X0: " << 0.5*h*x0 << std::endl;
					//std::cout << "X:  " << 0.5*h*x1 << std::endl;
				//}
				if (!success) {
					std::cout << "Fail with" << std::endl;
					std::cout << "M:  " << 0.5*h*mu1 << std::endl;
					std::cout << "X0: " << 0.5*h*x0 << std::endl;
					std::cout << "X:  " << 0.5*h*x1 << std::endl;
					return 0;
				}

				// Updating time velocity and momentum
				problem.mom_time(i,j,mu1);
				if (!success)
					x1 = x0;
				problem.vel_time(i,j,x1);

				// Updating position
				problem.pos(i+1,j,problem.pos(i,j)*Cay::eval(h*x1));
			}
		} // loop j

		problem.updateCSV(i,w_resample);
		displacement.insert(problem.pos(i,indice_ecoute).translationVector()[1]);
	} // loop i

	problem.endCSV();
	displacement.flush();
	displacement.close();
	
	//problem.writeCSVFile("results.csv",true,300);

	return 0;
}
