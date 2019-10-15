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

	Eigen::Matrix<double,6,1>
	dLdv0 (const Algebra g)
	{ return m_Inertia.asDiagonal()*g.toVector(); }

	Eigen::Matrix<double,6,1>
	dLdv1 (const Algebra g)
	{ return m_Constraint.asDiagonal()*(g.toVector()-Algebra::GeneratorVector(6)); }
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

	int n_time_steps = 20, n_space_steps = 11;
	problem.setSize(n_time_steps,n_space_steps);

	// linspace
	double h = 0.005, s = 0.1;
	problem.baselinstep(0.0,h,0.0,s);
	std::cout << "h: " << problem.step_size(0,0,0) << std::endl;
	std::cout << "s: " << problem.step_size(0,0,1) << std::endl;
	
	Eigen::Matrix<double,6,1> v_eta0;
	v_eta0 << 1.0, 1.5, 1.0, 0.0, 0.0, 1.0;
	Algebra eta0(v_eta0);

	Eigen::Matrix<double,6,1> v_eta1;
	v_eta1 << 1.004, 1.52, 1.005, -0.01, 0.0, 1.0;
	Algebra eta1(v_eta1);

	Group g0 = Group::Identity();
	problem.pos(0,0,g0);
	g0[5] = h;
	problem.pos(1,0,g0);

	int i,j;

	for (j=1; j<n_space_steps; j++) {
		problem.pos(0,j,problem.pos(0,j-1)*(s*eta0).cay());
		problem.pos(0,j,problem.pos(0,j-1)*(s*eta1).cay());
	}

	/* Call the integrator ***************************************************/

	BiVariational::Integrator<double,Group,Algebra> integrator(problem,comm);
	bool success = integrator.initialize();
	std::cout << "Init successful: " << success << std::endl;
	integrator.integrate();

	/* Exit ******************************************************************/

#ifdef HAVE_MPI
	MPI_Finalize();
#endif

	if (myRank == 0)
		std::cout << "End result: TEST PASSED" << std::endl;

	return 0;
}

