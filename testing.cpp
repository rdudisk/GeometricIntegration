#include <iostream>
#include <cmath>

#include "Geomi/Common"
#include "Geomi/Variational"


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

int
main (int argc, char* argv[])
{
	double h = 0.1;
	int n_steps = 25;

	RigidBody myProblem;
	myProblem.baselinstep(0.0,h,n_steps);

	Variational::Abstract::Integrator* integrator;
	Variational::CovariantStep<M,Group,Algebra>* step = new Variational::CovariantStep<M,Group,Algebra>(myProblem);
	integrator = new Variational::Integrator<M,Group,Variational::CovariantStepInternals<M,Group,Algebra>,Variational::Abstract::LieProblem<M,Group,Algebra>,Algebra>(myProblem, *step);

	Eigen::Matrix<double,3,1> Inertia(2.0/3.0,1.0,2.0);
	myProblem.Inertia(Inertia);

	Group pos0 = Group::Identity();
	myProblem.pos(0,pos0);

	Algebra vel0 = Algebra(cos(M_PI/3.0), 0.0, sin(M_PI/3.0));
	Group pos1 = step->posFromVel(h,pos0,vel0);
	myProblem.pos(1,pos1);

	integrator->initialize();
	integrator->integrate();
	
	myProblem.write2csv("results.csv");

	return 0;
}

