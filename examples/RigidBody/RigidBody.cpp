/**
 * This is an example of using the Runge-Kutta Munthe-Kaas module.
 *
 * The test problem is the free rigid body problem.
 * Reference can be found in ...
 *
 * \autor Pierre Carré
 */

#include <iostream>
#include <cmath>

#include "Geomi/Common"
#include "Geomi/RKMK"

#include <Eigen/Dense>

typedef SO3::Algebra<double> Algebra;

typedef NOXVector<3> Q;
typedef NOXVector<3> NOXV;

class RigidBodyProblem : public RKMK::Abstract::Problem<double,Q,Algebra>
{
	private:
		Eigen::Matrix<double,3,1> m_Iinv;
	
	public:
		RigidBodyProblem()
		{ m_Iinv << 1, 1, 1; }

		bool
		computeA (Algebra& A, const Q& x)
		{
			A = Algebra((-m_Iinv).asDiagonal()*x);
			return true;
		}

		bool
		computeJacobianA (std::vector<Algebra>& JA, const Q& x)
		{
			Eigen::Matrix<double,3,3> M = (-m_Iinv).asDiagonal();
			for (int i=0; i<3; i++)
				JA.push_back(Algebra(M.col(i)));
			return true;
		}
	
		Eigen::Matrix<double,3,1>&
		Iinv ()
		{ return m_Iinv; }

		void
		Iinv (Eigen::Matrix<double,3,1> val)
		{ m_Iinv = val; }
};

int
main (int argc, char* argv[])
{
	double h = 0.1;
	int n_steps = 1000;
	RigidBodyProblem myProblem;
	myProblem.baselinstep(0.0,h,n_steps);

	Eigen::Matrix<double,3,1> Iinv(3.0/2.0,1.0,0.5);
	myProblem.Iinv(Iinv);

	Eigen::Matrix<double,3,1> pos(cos(M_PI/3.0), 0.0, sin(M_PI/3.0));
	myProblem.pos(0,pos);

	RKMK::Abstract::Integrator& integrator = RKMK::Factory<double,Q,Algebra>::createIntegrator(myProblem,"RK 4");

	bool success = integrator.integrate();

	//myProblem.write2csv("results.csv");

	return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

