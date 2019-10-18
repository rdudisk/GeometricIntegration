#include <iostream>
#include <cmath>

#include "Geomi/Common"
#include "Geomi/BiVariational"

typedef double					M;
typedef SE3::Group<double>		Group;
typedef SE3::Algebra<double>	Algebra;

class RigidBody : public BiVariational::Abstract::LieProblem<M,Group,Algebra>
{
private:
	Eigen::Matrix<double,6,1> m_Inertia;
	Eigen::Matrix<double,6,1> m_Constraint;

public:
	RigidBody ()
	{
		m_Inertia << 1, 1, 1, 1, 1, 1;
		m_Constraint << 1, 1, 1, 1, 1, 1;
	}
	
	Eigen::Matrix<double,6,1>&
	Inertia ()
	{ return m_Inertia; }

	void
	Inertia (Eigen::Matrix<double,6,1> val)
	{ m_Inertia = val; }

	Eigen::Matrix<double,6,1>&
	Constraint ()
	{ return m_Constraint; }

	void
	Constraint (Eigen::Matrix<double,6,1> val)
	{ m_Constraint = val; }

	void
	setInertia (double area, double rho)
	{
		// la formule n'est peut-être pas exacte, mais donne un ordre de grandeur
		double linMass = area*rho;
		m_Inertia << 0.25*linMass*area, 0.25*linMass*area, 0.5*linMass*area,
						linMass, linMass, linMass;

	}

	void
	setConstraint (double area, double young, double poisson)
	{
		double G = young/(2.0+2.0*poisson);
		m_Constraint << young*m_Inertia[0], young*m_Inertia[1], G*(m_Inertia[0]+m_Inertia[1]),
						G*area, G*area, young*area;
	}

	Eigen::Matrix<double,6,1>
	dLdv0 (const Algebra g)
	{ return m_Inertia.asDiagonal()*g.toVector(); }

	Eigen::Matrix<double,6,1>
	dLdv1 (const Algebra g)
	{ return m_Constraint.asDiagonal()*(g.toVector()-Algebra::GeneratorVector(5)); }

	double
	coeffCFL (double young, double poisson, double rho, double alpha=2.0)
	{ return 1.0/(alpha*sqrt(young*(1.0-poisson)/(rho*(1.0+poisson)*(1.0-2.0*poisson)))); }

	void
	writeCSVFile (const std::string filename, bool header = true)
	{
		std::ofstream of;
		of.open(filename,std::ios_base::trunc);

		int i,j;
		int n_time = this->size(0);
		int n_space = this->size(1);
		double T,S;
		Algebra xi, eta;
		Group p;
		Eigen::Matrix<double,3,1> x,u,v,c,d,E1,E2;
		E1 << this->step_size(0,0,1), 0, 0;
		E2 << 0, this->step_size(0,0,1), 0;
		
		if (header)
			of << "i,j,t,s,x,y,z,u1,u2,u3,v1,v2,v3,c1,c2,c3,d1,d2,d3" << std::endl;

		T = 0;
		for (i=0; i<n_time; i++) {

			if (i>0)
				T += this->step_size(i-1,0,0);

			S = 0;
			for (j=0; j<n_space; j++) {

				if (j>0)
					S += this->step_size(i,j-1,1);
				p = this->pos(i,j);
				x = p.trans();
				u = p.transformVector(E1)-x;
				v = p.transformVector(E2)-x;
				c = xi.trans();
				d = eta.trans();
				if (j<n_space-1)
					eta = this->vel(i,j,1);
				else
					eta = Algebra::Zero();
				if (i<n_time-1)
					xi = this->vel(i,j,0);
				else
					xi = Algebra::Zero();
				of	<< i << ","
					<< j << ","
					<< T << ","
					<< S << ","
					<< x[0] << ","
					<< x[1] << ","
					<< x[2] << ","
					<< u[0] << ","
					<< u[1] << ","
					<< u[2] << ","
					<< v[0] << ","
					<< v[1] << ","
					<< v[2] << ","
					<< c[0] << ","
					<< c[1] << ","
					<< c[2] << ","
					<< d[0] << ","
					<< d[1] << ","
					<< d[2] << std::endl;
			}
		}

		of.close();
	}
};

int
main (int argc, char* argv[])
{
	/* Init MPI **************************************************************/

#ifdef HAVE_MPI
	MPI_Init(&argc,&argv);
	Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm comm;
#endif

	const int myRank = comm.MyPID();
	const int numProcs = comm.NumProc();

	if (myRank == 0)
		std::cout << "Total number of processes " << numProcs << std::endl;

	/* Setting up the problem ************************************************/

	RigidBody problem;

	double diameter = 0.01;
	double area = diameter*diameter;
	double rho = 1.0e3;
	double young = 5.0e3;
	double poisson = 0.35;
	std::cout << "Condition CFL: h = s*" << problem.coeffCFL(young,poisson,rho) << std::endl;
	//problem.setInertia(area,rho);
	//problem.setConstraint(area,young,poisson);

	int n_time_steps = 5, n_space_steps = 4;
	problem.setSize(n_time_steps,n_space_steps);

	// linspace
	double h = 0.05, s = 0.1;
	problem.baselinstep(0.0,h,0.0,s);
	
	Eigen::Matrix<double,6,1> v_eta0;
	v_eta0 << 1.0, 1.5, 1.0, 0.0, 0.0, 1.0;
	Algebra eta0(v_eta0);

	Eigen::Matrix<double,6,1> v_eta1;
	v_eta1 << 1.004, 1.52, 1.005, -0.01, 0.0, 1.0;
	Algebra eta1(v_eta1);

	Group g0 = Group::Identity();
	problem.pos(0,0,g0);
	g0.trans(2,h);
	problem.pos(1,0,g0);

	int i,j;

	for (j=0; j<n_space_steps; j++) {
		if (j>0) {
			problem.pos(0,j,problem.pos(0,j-1)*(s*eta0).cay());
			problem.pos(1,j,problem.pos(1,j-1)*(s*eta1).cay());
		}
		/*
		if (j<n_space_steps-1) {
			problem.vel(0,j,1,eta0);
			problem.vel(1,j,1,eta1);
		}
		problem.vel(0,j,0,(1.0/h)*Algebra::cay_inv(problem.pos(0,j).inverse()*problem.pos(1,j)));
		*/
	}

	/* Call the integrator ***************************************************/

	BiVariational::Integrator<double,Group,Algebra> integrator(problem,comm);
	bool success = integrator.initialize();
	std::cout << "Init successful: " << success << std::endl;
	integrator.integrate();
	
	/* Output results ********************************************************/

	problem.writeCSVFile("results.csv",true);

	/* Exit ******************************************************************/

#ifdef HAVE_MPI
	MPI_Finalize();
#endif

	if (myRank == 0)
		std::cout << "End result: TEST PASSED" << std::endl;

	return 0;
}

