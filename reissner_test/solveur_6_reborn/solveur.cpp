#include <iostream>
#include <cmath>
#include <string>

#include "header.hpp"

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

	double radius, rho, young, shear, k_factor, f_factor, pincement, tension, length_T0, h;
	int    n_space_steps, n_time_steps;

	pugi::xml_node xml_params = doc.child("config").child("params");
	young     = std::stod(xml_params.child("young").child_value());
	shear     = std::stod(xml_params.child("shear").child_value());
	rho       = std::stod(xml_params.child("rho").child_value());
	radius    = std::stod(xml_params.child("radius").child_value());
	k_factor  = std::stod(xml_params.child("kfactor").child_value());
	f_factor  = std::stod(xml_params.child("ffactor").child_value());
	pincement = std::stod(xml_params.child("pincement").child_value());
	tension   = std::stod(xml_params.child("tension").child_value());
	length_T0 = std::stod(xml_params.child("length").child_value());

	pugi::xml_node xml_integration = doc.child("config").child("integration");
	h = std::stod(xml_integration.child("time-step").child_value());
	/* n_space_steps = nombre d'intervalles
	 * le nombre de points spatiaux est donc n_space_steps+1
	 */
	n_space_steps = std::stoi(xml_integration.child("n-space-step").child_value());
	n_time_steps  = std::stoi(xml_integration.child("n-time-step").child_value());

	double l, length_default, area, alpha;
	area           = RigidBody::compute_area(radius);
	alpha          = tension/(young*area);
	length_default = length_T0/(1.0+alpha);
	l              = length_default/n_space_steps;
	//l              = length_T0/n_space_steps;
	std::cout << "L_t0: " << length_T0 << std::endl;
	std::cout << "L   : " << length_default << std::endl;
	std::cout << "Theoretical f0: " << (0.5/length_T0)*sqrt(tension/(rho*area*(1+alpha))) << std::endl;

	/* Checking up CFL condition requirement **********************************/

	double alpha_cfl = 0.1;
	double h_cfl=l*RigidBody::coeffCFL(young, shear, rho, alpha_cfl);
	
	if (h>=h_cfl) {
		std::string response;
		std::cout << "Warning: CFL condition not verified (factor=" << alpha_cfl << "):" << std::endl
			      << "\th = " << h << std::endl
				  << "\tlimit = " << h_cfl << std::endl;
		/*		  << "Do you want to continue ? [yn] ";
		std::cin >> response;
		if (response!="y")
			return 0;*/
	}

	/* Instantiate Reissner beam problem **************************************/

	RigidBody problem;

	problem.setTensors(radius,rho,young,shear);
	problem.setSize(n_space_steps+1);
	problem.setStepSize(h,l);
	problem.setCSV("results.csv");

	Displacement displacement_x(1.0/h,"displacement-x.csv");
	Displacement displacement_y(1.0/h,"displacement-y.csv");
	Displacement displacement_z(1.0/h,"displacement-z.csv");

	// TODO: clean that up
	double F_constant = f_factor/(rho*area);
	double F_duration = 0.01;
	Eigen::Matrix<double,6,1> Fvector;
	Fvector << 0.0, 0.0, 0.0, 0.0, F_constant, 0.0;
	//Fvector << 0.0, 0.0, 0.0, F_constant, 0.0, 0.0;
	int indice_pincement = (int) (pincement*((float)n_space_steps));
	int indice_ecoute = indice_pincement;
	Eigen::Matrix<double,3,1> listenVector;
	listenVector << 0.0, 0.0, 1.0;
	double elongation_factor = 1.0+(tension/(area*young));

	/* Inverse Legendre *******************************************************/

	InverseLegendre inv_leg(problem);
	RestrictedInverseLegendre restr_inv_leg(problem);

	/* Initialization *********************************************************/

	int i,j;
	bool verbose;
	bool success;
	double w_resample = (1.0/44100.0);
	
	Group   g0;
	Algebra x0,  x1,  e0,   e1;
	Vec6    mu0, mu1, sig0, sig1;
	Eigen::Matrix<double,6,1> E0, E4, eps_ref;

	// Setting up a equilibrium reference beam with null speed
	
	g0   = Group::Identity();
	E4   << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
	//E0   << 0.01, -0.03, 0.02, 0.98, 0.0, 0.01;
	//e0   = elongation_factor*Algebra(E4);
	//eps_ref = E4*length_default;
	eps_ref = E4;
	//e0   = Algebra(E0);
	x0   = Algebra::Zero();
	x1   = Algebra::Zero();
	//x1[4]= 0.01*h;
	mu0  = Vec6::Zero();
	sig0 = Vec6::Zero();

	// init position
	problem.pos(0,0,g0);
	for (j=1; j<=n_space_steps; j++) {
		g0.translation(0,length_T0*float(j)/float(n_space_steps));
		problem.pos(0,j,g0);
		//problem.pos(0,j,problem.pos(0,j-1)*Cay::eval(l*e0));
	}

	// init time velocity and momentum
	for (j=0; j<=n_space_steps; j++) {
		if (false) { //(j==3) {
			problem.vel_time(0,j,x1);
			problem.mom_time(0,j,Cay::inv_right_diff_star(h*x1,Algebra(problem.Inertia()*(x1.vector()))).vector());
		}
		else {
			problem.vel_time(0,j,x0);
			problem.mom_time(0,j,mu0);
		}
	}

	// init space velocity and momentum
	for (j=0; j<n_space_steps; j++) {
		e0   = (1.0/l)*Cay::inv(problem.pos(0,j).inverse()*problem.pos(0,j+1));
		sig0 = (-1.0)*Cay::inv_right_diff_star(l*e0,Algebra(problem.Constraint()*(e0.vector()-eps_ref))).vector();
		problem.vel_space(0,j,e0);
		problem.mom_space(0,j,sig0);
	}
	
	double displacement_x0 = problem.pos(0,indice_ecoute).translationVector()[0];
	displacement_x.insert(problem.pos(0,indice_ecoute).translationVector()[0]-displacement_x0);
	displacement_y.insert(problem.pos(0,indice_ecoute).translationVector()[1]);
	displacement_z.insert(problem.pos(0,indice_ecoute).translationVector()[2]);
	//displacement.insert(1-problem.pos(0,indice_ecoute).rotateVector(listenVector)[2]);


	// Updating position for i=1
	for (j=0; j<=n_space_steps; j++) {
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
		for (j=0; j<=n_space_steps; j++) {
			if (j<n_space_steps) {
				// Updating space velocity and momentum
				e1   = (1.0/l)*(Cay::inv(problem.pos(i,j).inverse()*problem.pos(i,j+1)));
				sig1 = (-1.0)*(Cay::inv_right_diff_star(l*e1,Algebra(problem.Constraint()*(e1.vector()-eps_ref))).vector());
				problem.vel_space(i,j,e1);
				problem.mom_space(i,j,sig1);
			}
			//else
				//e1 = Algebra(E4);

			if (j==0) {
				x0   = problem.vel_time(i-1,j);
				mu0  = problem.mom_time(i-1,j);
				sig1 = problem.mom_space(i,j);

				mu1 =  Algebra::static_Ad_star(Cay::eval(h*x0),Algebra(mu0)).vector();
				mu1 += (-2.0*h/l)*sig1;

				// solve for x1
				restr_inv_leg.setData(h,mu1);
				success = restr_inv_leg.computeSolution();
				//inv_leg.setData(h,mu1);
				//success = inv_leg.computeSolution();

				x1  = restr_inv_leg.getAlgebraSolution();
				mu1 = (Cay::inv_right_diff_star(h*x1,Algebra(problem.Inertia()*x1.vector()))).vector();
				//x1  = inv_leg.getAlgebraSolution();

				// Updating time velocity and momentum
				problem.mom_time(i,j,mu1);
				if (!success) {
					x1 = x0;
					return 0;
				}	
				problem.vel_time(i,j,x1);

				// Updating position
				problem.pos(i+1,j,problem.pos(i,j)*Cay::eval(h*x1));
			}
			else if (j==n_space_steps) {
				x0   = problem.vel_time(i-1,j);
				mu0  = problem.mom_time(i-1,j);
				e0   = problem.vel_space(i,j-1);
				sig0 = problem.mom_space(i,j-1);

				mu1 =  Algebra::static_Ad_star(Cay::eval(h*x0),Algebra(mu0)).vector();
				mu1 += (2.0*h/l)*(Algebra::static_Ad_star(Cay::eval(l*e0),Algebra(sig0)).vector());

				//mu1[0] = 0.0;
				//mu1[3] = 0.0; mu1[4] = 0.0; mu1[5] = 0.0;

				// solve for x1
				restr_inv_leg.setData(h,mu1);
				success = restr_inv_leg.computeSolution();
				//inv_leg.setData(h,mu1);
				//success = inv_leg.computeSolution();

				x1  = restr_inv_leg.getAlgebraSolution();
				mu1 = (Cay::inv_right_diff_star(h*x1,Algebra(problem.Inertia()*x1.vector()))).vector();
				//x1  = inv_leg.getAlgebraSolution();

				// Updating time velocity and momentum
				problem.mom_time(i,j,mu1);
				if (!success) {
					x1 = x0;
					return 0;
				}	
				problem.vel_time(i,j,x1);

				// Updating position
				problem.pos(i+1,j,problem.pos(i,j)*Cay::eval(h*x1));
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
					mu1 += h*(i_double*h/F_duration)*Fvector;

				inv_leg.setData(h,mu1);
				success = inv_leg.computeSolution();
				x1 = inv_leg.getAlgebraSolution();

				// Updating time velocity and momentum
				problem.mom_time(i,j,mu1);
				if (!success) {
					x1 = x0;
					return 0;
				}	
				problem.vel_time(i,j,x1);

				// Updating position
				problem.pos(i+1,j,problem.pos(i,j)*Cay::eval(h*x1));
			}
		} // loop j

		problem.updateCSV(i,w_resample);
		displacement_x.insert(problem.pos(i,indice_ecoute).translationVector()[0]-displacement_x0);
		displacement_y.insert(problem.pos(i,indice_ecoute).translationVector()[1]);
		displacement_z.insert(problem.pos(i,indice_ecoute).translationVector()[2]);
		//displacement.insert(1-problem.pos(0,indice_ecoute).rotateVector(listenVector)[2]);
	} // loop i

	problem.endCSV();
	displacement_x.flush();
	displacement_x.close();
	displacement_y.flush();
	displacement_y.close();
	displacement_z.flush();
	displacement_z.close();
	
	return 0;
}
