#ifndef CLASSES_HPP
#define CLASSES_HPP

#include "Geomi/src/Common/utils.hpp"
#include "Geomi/src/Common/NOXVector.hpp"
#include "Geomi/src/Common/Abstract_NOXStep.hpp"
#include "Geomi/src/Common/NOXGroup.hpp"
#include "Geomi/src/Common/LieGroupBase.hpp"
#include "Geomi/src/Common/LieAlgebraBase.hpp"
#include "Geomi/src/Common/SO3_Group.hpp"
#include "Geomi/src/Common/SO3_Algebra.hpp"
#include "Geomi/src/Common/SE3_Group.hpp"
#include "Geomi/src/Common/SE3_Algebra.hpp"
//#include "src/Common/csv.hpp"
//#include "src/Common/Syst.hpp"
//#include "src/Common/DiscSyst.hpp"
//#include "src/Common/GaussLegendre.hpp"
//#include "src/Common/LagrangeInterpolation.hpp"

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

typedef double					M;
typedef SE3::Group<double>		Group;
typedef SE3::Algebra<double>	Algebra;
typedef Eigen::Matrix<double,6,1>	Vec6;
typedef Eigen::Matrix<double,3,1>	Vec3;

class MultiSyst;
class DiscMultiSyst;
class RigidBody;
class SolveMe;

class MultiSyst
{
protected:
	M		m_base_time;
	M		m_base_space;
	Group	m_pos;
	Algebra m_vel_time;
	Algebra m_vel_space;
	Vec6	m_mom_time;
	Vec6	m_mom_space;

public:
	MultiSyst ( ) { }
	MultiSyst (M base_time, M base_space) : m_base_time(base_time), m_base_space(base_space) { }
	~MultiSyst ( ) { }
	M base_time (void) const;
	void base_time(M b);
	M base_space (void) const;
	void base_space(M b);
	Group pos ( ) const;
	void pos (Group p);
	Algebra vel_time ( ) const;
	void vel_time (Algebra v);
	Algebra vel_space ( ) const;
	void vel_space (Algebra v);
	Vec6 mom_time ( ) const;
	void mom_time (Vec6 m);
	Vec6 mom_space ( ) const;
	void mom_space (Vec6 m);
};


class DiscMultiSyst
{
protected:
	std::vector<MultiSyst> m_node;
	size_t m_size[2];

public:
	DiscMultiSyst ( )
	{
		m_size[0] = 0;
		m_size[1] = 0;
		m_node.clear();
	}
	~DiscMultiSyst ( ) { }

	void setSize (const size_t i, const size_t j);
	MultiSyst const& operator[] (size_t index) const;
	MultiSyst& operator[] (size_t index);
	size_t size(const size_t dir) const;
	const size_t getIndex (const size_t i, const size_t j) const;
	M step_size (const size_t dir) const;
	M base_time (const size_t& i) const;
	M base_space (const size_t& j) const;
	Group pos (const size_t& i, const size_t& j) const;
	void pos (const size_t& i, const size_t& j, const Group& p);
	Algebra vel_time (const size_t& i, const size_t& j) const;
	void vel_time (const size_t& i, const size_t& j, const Algebra& v);
	Algebra vel_space (const size_t& i, const size_t& j) const;
	void vel_space (const size_t& i, const size_t& j, const Algebra& v);
	Vec6 mom_time (const size_t& i, const size_t& j) const;
	void mom_time (const size_t& i, const size_t& j, const Vec6& m);
	Vec6 mom_space (const size_t& i, const size_t& j) const;
	void mom_space (const size_t& i, const size_t& j, const Vec6& m);
	static const unsigned int dof ( );
	void baselinstep (M t_inf_lim, M t_step_size, M s_inf_lim, M s_step_size);
};

class RigidBody : public DiscMultiSyst
{
private:
	Eigen::Matrix<double,6,6> m_Inertia;
	Eigen::Matrix<double,6,6> m_Constraint;

public:
	RigidBody ()
	{ m_Inertia = Eigen::Matrix<double,6,6>::Identity; m_Constraint = m_Inertia; }
	
	Eigen::Matrix<double,6,6>& Inertia ();
	void Inertia (Eigen::Matrix<double,6,6> val);
	Eigen::Matrix<double,6,6>& Constraint ();
	void Constraint (Eigen::Matrix<double,6,6> val);
	void setInertia (double area, double rho);
	void setConstraint (double area, double young, double poisson);
	double coeffCFL (double young, double poisson, double rho, double alpha=3.0);
	void writeCSVFile (const std::string filename, bool header = true);
};

class SolveMe : public ::Abstract::NOXStep<NOXVector<6>,1>
{
protected:
	double m_h;
	Vec6 m_mu;
	Vec6 m_x0;

	RigidBody& m_problem;

public:
	SolveMe (RigidBody& problem)
	: m_problem(problem) { }

	void setData (double h, Vec6 mu);
	void setx0 (Vec6 x0);
	double h () const;
	Vec6 mu () const;
	const NOXVector<6> getInitialGuess ();
	bool computeF (NOXVector<6>& f, const NOXVector<6>& x);
	bool computeJacobian (Eigen::Matrix<double,6,6>& J, const NOXVector<6>& x);
};

#endif
