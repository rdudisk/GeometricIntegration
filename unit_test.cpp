#define BOOST_TEST_MODULE "SO(3) group test"
#include <boost/test/included/unit_test.hpp>

#include "headers/SO3Group.hpp"
#include "headers/SO3Algebra.hpp"

BOOST_AUTO_TEST_CASE ( so3_basics )
{
	using namespace Lie::SO3;
	typedef Group<float> Group_f;
	typedef Algebra<float> Algebra_f;

	Group_f Id_g = Group_f::Identity();
	BOOST_TEST( Id_g.isApprox(Group_f(Id_g.q())) );

	Algebra_f Id_a = Algebra_f::Identity();
	BOOST_TEST( Id_g.isApprox(Id_a.exp()) );
}
