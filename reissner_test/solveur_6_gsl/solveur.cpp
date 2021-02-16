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

	double diameter, rho, young, poisson, area, k_factor, f_factor;

	pugi::xml_node xml_params = doc.child("config").child("params");
	young    = std::stod(xml_params.child("young").child_value());
	poisson  = std::stod(xml_params.child("poisson").child_value());
	rho      = std::stod(xml_params.child("rho").child_value());
	diameter = std::stod(xml_params.child("diameter").child_value());
	area     = diameter*diameter;
	k_factor = std::stod(xml_params.child("kfactor").child_value());
	f_factor = std::stod(xml_params.child("ffactor").child_value());

	double l, h;
	int n_space_steps, n_time_steps;

	pugi::xml_node xml_integration = doc.child("config").child("integration");
	l = std::stod(xml_integration.child("space-step").child_value());
	h = std::stod(xml_integration.child("time-step").child_value());
	n_space_steps = std::stoi(xml_integration.child("n-space-step").child_value());
	n_time_steps  = std::stoi(xml_integration.child("n-time-step").child_value());

	/* Checking up CFL condition requirement **********************************/

	double cfl, cfl_factor;
	double lame_lambda, lame_mu;
	lame_lambda = young*poisson/((1.0+poisson)*(1.0-2.0*poisson));
	lame_mu     = young/(2.0*(1.0+poisson));
	cfl_factor  = 0.1;
	cfl	        = cfl_factor*l/sqrt((lame_lambda+2.0*lame_mu)/rho);
	
	if (h>=cfl) {
		std::string response;
		std::cout << "Warning: CFL condition not verified (factor=0.1):" << std::endl
			      << "\th = " << h << std::endl
				  << "\tlimit = " << cfl << std::endl
				  << "Do you want to continue ? [yn] ";
		std::cin >> response;
		if (response!="y")
			return 0;
	}

	/* Instantiate Reissner beam problem **************************************/

	RigidBody problem;

	problem.setInertia(area,rho);
	problem.setConstraint(area,young,poisson);
	problem.setSize(n_time_steps,n_space_steps);
	problem.baselinstep(0.0,h,0.0,l);

	/* Initialization *********************************************************/

	int i,j;
	bool verbose;
	int success;
	
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
		if (j==3) {
			problem.vel_time(0,j,x1);
			problem.mom_time(0,j,Cay::inv_right_diff_star(h*x1,Algebra(problem.Inertia()*(x1.vector()))).vector());
		}
		else {
			problem.vel_time(0,j,x0);
			problem.mom_time(0,j,mu0);
		}
		if (j<n_space_steps-1) {
			problem.vel_space(0,j,e0);
			problem.mom_space(0,j,sig0);
		}
	}

	// Updating position for i=1
	
	for (j=0; j<n_space_steps; j++) {
		problem.pos(1,j,problem.pos(0,j)*Cay::eval(h*problem.vel_time(0,j)));
	}

	/* Integration loop *******************************************************/

	// TODO: clean that up
	
	int solve_err_nb = 0, max_err_nb = 20;
	Algebra* x_sol = new Algebra;
	Eigen::Matrix<double,6,6> Kmatrix;
	Eigen::Matrix<double,6,1> Kvector;
	Kvector << 0.0, 0.0, 0.0, 1.0, 1.0, 1.0;
	Kmatrix = k_factor*Kvector.asDiagonal();
	double F_constant = f_factor/(rho*area);
	double F_duration = 0.01;
	Eigen::Matrix<double,6,1> Fvector;
	Fvector << 0.0, 0.0, 0.0, 0.0, F_constant, 0.0;
	double i_double;
	
	for (i=1; i<n_time_steps-1; i++) {
		for (j=0; j<n_space_steps; j++) {
			if (j<n_space_steps-1) {
				// Updating space velocity and momentum
				e1   = (1.0/l)*Cay::inv(problem.pos(i,j).inverse()*problem.pos(i,j+1));
				sig1 = (-1.0)*Cay::inv_right_diff_star(l*e1,Algebra(problem.Constraint()*(e1.vector()-E4))).vector();
				problem.vel_space(i,j,e1);
				problem.mom_space(i,j,sig1);
			}
			//else
				//e1 = Algebra(E4);

			if (j==0 || j==(n_space_steps-1)) {
				// Updating time velocity and momentum
				x1  = Algebra::Zero();
				mu1 = Vec6::Zero();
				problem.vel_time(i,j,x1);
				problem.mom_time(i,j,mu1);

				// Updating position
				problem.pos(i+1,j,problem.pos(i,j));
			}
			else {
				x0   = problem.vel_time(i-1,j);
				mu0  = problem.mom_time(i-1,j);
				e0   = problem.vel_space(i,j-1);
				sig0 = problem.mom_space(i,j-1);

				// Writing DEL equation
				mu1 = Algebra::static_Ad_star(Cay::eval(h*x0),Algebra(mu0)).vector();
				mu1 -= (h/l)*(Algebra::static_Ad_star(Cay::eval(l*e0),Algebra(sig0)).vector() -sig1);
				// Adding friction forces
				mu1 -= h*l*Kmatrix*x0.vector();
				// Adding external control forces
				i_double = (double) i;
				//if (j==3 && i_double*h < F_duration)
					//mu1 += h*l*(i_double*h/F_duration)*Fvector;

				// Solving for x1
				std::cout << "Solving node " << i << ", " << j << "... ";
				success = solve_speed(mu1,h,x0,problem.Inertia(),x_sol);
				solve_err_nb += success;
				if (!success) std::cout << "success!" << std::endl;
				else std::cout << "fail!" << std::endl;
				x1 = *x_sol;

				// Updating time velocity and momentum
				problem.vel_time(i,j,x1);
				problem.mom_time(i,j,mu1);

				// Updating position
				problem.pos(i+1,j,problem.pos(i,j)*Cay::eval(h*x1));
			}
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

	problem.writeCSVFile("results.csv",true,300);

	//std::cout << "End result: TEST PASSED" << std::endl;

	return 0;
}

	/* Garder: ça marche bien pour les parametres donnes par Demoures
	Eigen::Matrix<double,6,1> v_e0, v_e1, E4;
	v_e0 << 1.0, 1.0, 1.5, 1.0, 0.0, 0.0;
	v_e1 << 1.005, 1.004, 1.52, 1.0, 0.0, -0.01;
	E4   << 0, 0, 0, 1, 0, 0;
	e0   = Algebra(v_e0);
	e1   = Algebra(v_e1);

	g0 = Group::Identity();
	problem.pos(0,0,g0);
	g0.translation(2,h);
	problem.pos(1,0,g0);


	for (j=0; j<n_space_steps; j++) {
		if (j>0) {
			problem.pos(0,j,problem.pos(0,j-1)*Cay::eval(l*e0));
			problem.pos(1,j,problem.pos(1,j-1)*Cay::eval(l*e1));
		}
		x0  = (1.0/h)*Cay::inv(problem.pos(0,j).inverse()*problem.pos(1,j));
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
	*/

