#include <iostream>
#include <cmath>

#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/ellint_1.hpp>

#include "Geomi/Common"
#include "Geomi/Variational"


#define N_BODIES				2
#define BODY_DOF				3
typedef double					M;
typedef NOXVector<BODY_DOF*N_BODIES>	Q;
typedef NOXVector<BODY_DOF*N_BODIES>	TQ;

/*
 * Probleme tres simple : L(q,v)=1/2*m*v^2
 * Solution : q,v = cste
 */
typedef NOXVector<2>	N2;

class MyProblem : public Variational::Abstract::Problem<M,N2>
{
private:

public:
	Eigen::Matrix<double,N2::DOF,1>
	dLdq (const N2 q, const N2 v)
	{
		return Eigen::Matrix<double,N2::DOF,1>::Zero();
	}

	Eigen::Matrix<double,N2::DOF,1>
	dLdv (const N2 q, const N2 v)
	{
		return Eigen::Matrix<double,N2::DOF,1>(v);
	}

	Eigen::Matrix<double,N2::DOF,N2::DOF>
	JqdLdq (const N2 q, const N2 v)
	{
		return Eigen::Matrix<double,N2::DOF,N2::DOF>::Zero();
	}

	Eigen::Matrix<double,N2::DOF,N2::DOF>
	JvdLdq (const N2 q, const N2 v)
	{
		return Eigen::Matrix<double,N2::DOF,N2::DOF>::Zero();
	}

	Eigen::Matrix<double,N2::DOF,N2::DOF>
	JqdLdv (const N2 q, const N2 v)
	{
		return Eigen::Matrix<double,N2::DOF,N2::DOF>::Zero();
	}

	Eigen::Matrix<double,N2::DOF,N2::DOF>
	JvdLdv (const N2 q, const N2 v)
	{
		return Eigen::Matrix<double,N2::DOF,N2::DOF>::Identity();
	}
};

/**
 * 2 bodies problem
 * L(q,v) = (1/2)*(m1*v1**2+m2*v2**2) + G*m1*m2/||q1-q2||**2
 * where q = (q1,q2)^T and v = (v1,v2)^T
 */
class KeplerProblem : public Variational::Abstract::Problem<M,Q>
{
private:
	double m_G;
	double m_m1;
	double m_m2;

public:
	void
	G (double _G)
	{
		m_G = _G;
	}

	void
	m1 (double _m1)
	{
		m_m1 = _m1;
	}

	void
	m2 (double _m2)
	{
		m_m2 = _m2;
	}

	Eigen::Matrix<double,Q::DOF,1>
	dLdq (const Q q, const Q v)
	{
		Eigen::Matrix<double,BODY_DOF,1> q1 = q.head<3>(), q2 = q.tail<3>();

		Eigen::Matrix<double,Q::DOF,1> _dLdq;

		_dLdq.head<3>() = -2.0*m_G*m_m1*m_m2*pow((q1-q2).norm(),-4)*(q1-q2);
		_dLdq.tail<3>() = -_dLdq.head<3>();

		return _dLdq;
	}

	Eigen::Matrix<double,Q::DOF,1>
	dLdv (const Q q, const Q v)
	{
		Eigen::Matrix<double,BODY_DOF,1> v1 = v.head<3>(), v2 = v.tail<3>();

		Eigen::Matrix<double,Q::DOF,1> _dLdv;

		_dLdv.head<3>() = m_m1*v1;	
		_dLdv.tail<3>() = m_m2*v2;	

		return _dLdv;
	}

	Eigen::Matrix<double,Q::DOF,Q::DOF>
	JqdLdq (const Q q, const Q v)
	{
		Eigen::Matrix<double,BODY_DOF,1> q1 = q.head<3>(), q2 = q.tail<3>();

		Eigen::Matrix<double,Q::DOF,Q::DOF> _JqdLdq;

		double invSquareNorm = 1.0/((q1-q2).dot(q1-q2));
		double frontConst = m_G*m_m1*m_m2;
		
		_JqdLdq.block<3,3>(0,0) = -2.0*frontConst*pow(invSquareNorm,2)*Eigen::Matrix<double,3,3>::Identity()
									+ 8.0*frontConst*pow(invSquareNorm,3)*((q1-q2)*(q1-q2).transpose());
		_JqdLdq.block<3,3>(3,3) = _JqdLdq.block<3,3>(0,0);
		_JqdLdq.block<3,3>(3,0) = -_JqdLdq.block<3,3>(0,0);
		_JqdLdq.block<3,3>(0,3) = -_JqdLdq.block<3,3>(0,0);

		return _JqdLdq;
	}

	Eigen::Matrix<double,Q::DOF,Q::DOF>
	JvdLdq (const Q q, const Q v)
	{
		return Eigen::Matrix<double,Q::DOF,Q::DOF>::Zero();
	}

	Eigen::Matrix<double,Q::DOF,Q::DOF>
	JqdLdv (const Q q, const Q v)
	{
		return Eigen::Matrix<double,Q::DOF,Q::DOF>::Zero();
	}

	Eigen::Matrix<double,Q::DOF,Q::DOF>
	JvdLdv (const Q q, const Q v)
	{
		Eigen::Matrix<double,Q::DOF,Q::DOF> _JvdLdv;
		
		_JvdLdv.block<3,3>(0,0) = m_m1*Eigen::Matrix<double,3,3>::Identity();
		_JvdLdv.block<3,3>(3,3) = m_m2*Eigen::Matrix<double,3,3>::Identity();
		_JvdLdv.block<3,3>(0,3) = Eigen::Matrix<double,3,3>::Zero();
		_JvdLdv.block<3,3>(3,0) = Eigen::Matrix<double,3,3>::Zero();

		return _JvdLdv;
	}

};

int
main (int argc, char* argv[])
{
	double h = 0.1;
	int n_steps = 4;

	/*
	 * Initial data:
	 *				| Sun					| Jupiter				|
	 * mass			| 1.00000597682			| 0.0000954786104043	|
	 * q0(1)		| 0.0					| -3.5023653			|
	 * q0(2)		| 0.0					| -3.8169847			|
	 * q0(3)		| 0.0					| -1.5507963			|
	 * v0(1)		| 0.0					| 0.00565429			|
	 * v0(2)		| 0.0					| -0.00412490			|
	 * v0(3)		| 0.0					| -0.00190589			|
	 */
	KeplerProblem myProblem;
	myProblem.baselinstep(0.0,h,n_steps);
	Q pos, vel;
	myProblem.G(2.95912208286e-4);
	myProblem.m1(1.00000597682);
	myProblem.m2(0.0000954786104043);
	pos << 0.0, 0.0, 0.0, -3.5023653, -3.8169847, -1.5507963;
	vel << 0.0, 0.0, 0.0, 0.00565429, -0.00412490, -0.00190589;
	myProblem.pos(0,pos);

	// for now this is the only way to initialize the 2nd position consistently with the step definition
	// TODO !!
	myProblem.pos(1,pos+h*vel);

	//Variational::Abstract::Integrator& integrator = Variational::Factory<double,Q,TQ>::createIntegrator(myProblem,"Midpoint");
	
	Variational::Abstract::Integrator& integrator = Variational::Factory<double,Q,TQ>::createIntegrator(myProblem,"Galerkin P2N2Gau");


	/*
	MyProblem myProblem;
	myProblem.baselinstep(0.0,h,n_steps);
	N2 pos, vel;
	pos << 0, 0;
	vel << 1, -1;
	myProblem.pos(0,pos);
	myProblem.pos(1,pos+h*vel);

	Variational::Abstract::Integrator& integrator = Variational::Factory<double,N2,N2>::createIntegrator(myProblem,"Explicit Euler");
	*/

	integrator.integrate();
	myProblem.write2csv("res_midpoint2.csv");

	return 0;
}

