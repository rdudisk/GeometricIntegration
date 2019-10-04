#include <iostream>
#include <cmath>

#include "Geomi/Common"
#include "Geomi/MultiVariational"

typedef double					M;
typedef SO3::Group<double>		Group;
typedef SO3::Algebra<double>	Algebra;

class RigidBody : public MultiVariational::Abstract::LieProblem<Group,Algebra,2>
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

	/*
	template <>
	Eigen::Matrix<double,3,1>
	dLdv<0> (const Algebra g)
	{ return m_Inertia.asDiagonal()*g.toVector(); }

	template <>
	Eigen::Matrix<double,3,1>
	dLdv<1> (const Algebra g)
	{ return m_Inertia.asDiagonal()*g.toVector(); }

	Eigen::Matrix<double,3,3>
	JvdLdv<0,0> (const Algebra)
	{ return m_Inertia.asDiagonal(); }

	Eigen::Matrix<double,3,3>
	JvdLdv<0,1> (const Algebra)
	{ return m_Inertia.asDiagonal(); }

	Eigen::Matrix<double,3,3>
	JvdLdv<1,0> (const Algebra)
	{ return m_Inertia.asDiagonal(); }

	Eigen::Matrix<double,3,3>
	JvdLdv<1,1> (const Algebra)
	{ return m_Inertia.asDiagonal(); }
	*/
};
/*
template <>
Eigen::Matrix<double,3,1>
RigidBody::dLdv<0> (const Algebra g)
{ return m_Inertia.asDiagonal()*g.toVector(); }
*/

namespace ns {
	template <typename U>
	class A {
		public:
		template <typename V>
		void foo(V);
	};

	class B : public A<int>
	{
	};

	//template <>
	//template <>
	//void A<int>::foo<double> (double);
	
	//void
	//B::template foo<double>(double);
} // namespace

//template <>
//template <>
//void ns::A<int>::foo<double> (double d) { std::cout << d << std::endl; };

//void ns::B::template foo<double> (double d) { std::cout << d << std::endl; };

int
main (int argc, char* argv[])
{
	return 0;
}

