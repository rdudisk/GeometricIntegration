#include <iostream>
#include <cmath>

//#include "Geomi/Common"
//#include "Geomi/BiVariational"

#include <Epetra_config.h>

#ifdef HAVE_MPI
	#include <mpi.h>
	#include <Epetra_MpiComm.h>
#else
	#include <Epetra_SerialComm.h>
#endif


/*
typedef double					M;
typedef SO3::Group<double>		Group;
typedef SO3::Algebra<double>	Algebra;

class RigidBody : public Variational::Abstract::LieProblem<M,Group,Algebra>
{
private:
	Eigen::Matrix<double,3,1> m_Inertia;

public:
	RigidBody ()
	{ m_Inertia << 1, 1, 1; }
	
	Eigen::Matrix<double,3,1>&
	Inertia ()
	{ return m_Inertia; }

	void
	Inertia (Eigen::Matrix<double,3,1> val)
	{ m_Inertia = val; }

	Eigen::Matrix<double,3,1>
	dLdv (const Algebra g)
	{ return m_Inertia.asDiagonal()*g.toVector(); }

	Eigen::Matrix<double,3,3>
	JvdLdv (const Algebra)
	{ return m_Inertia.asDiagonal(); }
};
*/

void
routine (const Epetra_Comm &comm, std::ostream &out)
{
	out << "MyPID = " << comm.MyPID() << std::endl;
}

int
main (int argc, char* argv[])
{

#ifdef HAVE_MPI
	MPI_Init(&argc,&argv);
	Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm comm;
#endif

	const int myRank = comm.MyPID();
	const int numProcs = comm.NumProc();

	if (myRank > 0)
		std::cout << "Total number of processes " << numProcs << std::endl;

	routine(comm,std::cout);

	if (comm.MyPID() == 0)
		std::cout << "End result: TEST PASSED" << std::endl;

#ifdef HAVE_MPI
	MPI_Finalize();
#endif

	return 0;
}

