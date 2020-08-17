#define BOOST_TEST_MODULE "Geomi Common module test"
#include <boost/test/included/unit_test.hpp>

#include <vector>

#include "Geomi/Common"

#include "lie/lie_SE3_operator.hpp"

BOOST_AUTO_TEST_CASE ( common, * boost::unit_test::tolerance(1e-12) )
{
	int i,j;

	/* Lagrange polynomials tests */

	LagrangePolynomials lp;

	std::vector<double> dates;
	dates.push_back(0.0);
	dates.push_back(0.5);
	dates.push_back(1.0);

	std::vector<std::vector<double>> poly = lp.polynomials(0,dates);
	std::vector<std::vector<double>> poly_der = lp.polynomials_derivatives(0,dates);

	BOOST_TEST ( poly.size()==dates.size() );
	BOOST_TEST ( poly[0].size()==1 );
	BOOST_TEST ( poly[0][0]==1.0 );
	BOOST_TEST ( poly[2][0]==1.0 );

	BOOST_TEST ( poly_der.size()==dates.size() );
	BOOST_TEST ( poly_der[0].size()==1 );
	BOOST_TEST ( poly_der[0][0]==0.0 );
	BOOST_TEST ( poly_der[2][0]==0.0 );

	poly = lp.polynomials(2,dates);

	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			if (i==j)
				BOOST_TEST ( poly[i][j]==1.0 );
			else
				BOOST_TEST ( poly[i][j]==0.0 );
		}
	}

	/* SO(3) tests */

	typedef SO3::Group<double> SO3G;
	typedef SO3::Algebra<double> SO3A;

	SO3A w1(1.1,2.8,1.9);
	SO3A w2 = SO3A::cay_inv(w1.cay());

	for (i=0; i<3; i++)
		BOOST_TEST ( w1[i] == w2[i] );

	/* SE(3) tests ***********************************************************/

	typedef SE3::Group<double> SE3G;
	typedef SE3::Algebra<double> SE3A;

	typedef Lie::SE3::Group<double> SE3Gok;
	typedef Lie::SE3::Algebra<double> SE3Aok;
	typedef Lie::SE3::Operator<double> SE3Ook;

	SE3G g1, g2;
	SE3A a1, a2, a3, a4;
	SE3Gok G1, G2;
	SE3Aok A1, A2, A3, A4;

	/* Group inverse */
	g1 = SE3G::Random();
	g2 = g1*(g1.inverse());

	BOOST_TEST ( g2.isApprox(SE3G::Identity()) );

	/* 3D vector transformation */

	Eigen::Matrix<double,3,1> v3D1 = Eigen::Matrix<double,3,1>::Random().normalized();
	Eigen::Matrix<double,3,1> v3D2 = g1.transformVector(v3D1) - g1.trans();
	Eigen::Matrix<double,3,1> v3D3 = g1.rotateVector(v3D1);

	BOOST_TEST ( v3D1.norm() == v3D2.norm() );
	BOOST_TEST ( v3D1.norm() == v3D3.norm() );

	for (i=0; i<3; i++)
		BOOST_TEST ( v3D2[i] == v3D3[i] );

	/* Algebra inverse */
	a1 = SE3A::Random();
	a2 = a1+(a1.inverse());

	for (i=0; i<6; i++)
		BOOST_TEST ( a2[i] == 0.0 );

	/* Ad_star */
	/* Test the property Ad^*_{g^{-1}}\circ Ad^*_g = Ad^*_e */
	a1 = SE3A::Random();
	g1 = SE3G::Random();
	g2 = g1.inverse();
	a3 = SE3A::static_Ad_star(g1,SE3A::static_Ad_star(g2,a1));
	a4 = SE3A::static_Ad_star(SE3G::Identity(),a1);

	for (i=0; i<6; i++)
		BOOST_TEST ( a3[i] == a4[i] );

	/* Cayley and inverse */
	a1 = SE3A::Random();
	a2 = SE3A::cay_inv(a1.cay());

	for (i=0; i<6; i++)
		BOOST_TEST ( a1[i] == a2[i] );

	g1 = SE3A::Zero().cay();

	for (i=0; i<6; i++)
		BOOST_TEST ( g1.isApprox(SE3G::Identity()) );

	/* Cayley */

	Eigen::Matrix<double,6,1> v = Eigen::Matrix<double,6,1>::Random().normalized();
	a1 = SE3A(v);
	A1 = SE3Aok(v);

	for (i=0; i<6; i++)
		BOOST_TEST ( a1[i] == A1.vect()[i] );

	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			BOOST_TEST ( (a1.toMatrix())(i,j) == (A1.matrix())(i,j) );

	g1 = a1.cay();
	G1 = SE3Ook::cay(A1);
	
	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			BOOST_TEST ( (g1.matrix())(i,j) == (G1.matrix())(i,j) );

	a1 = SE3A::cay_inv(g1);
	A1 = SE3Ook::cay_inv(G1);

	for (i=0; i<6; i++)
		BOOST_TEST ( a1[i] == A1.vect()[i] );

	a1 = SE3A(v);
	A1 = SE3Aok(v);
	v = Eigen::Matrix<double,6,1>::Random().normalized();
	a2 = SE3A(v);
	A2 = SE3Aok(v);

	a3 = a1.dCayRInvStar(a2);
	A3 = SE3Ook::dcay_inv_star(A1,A2);

	for (i=0; i<6; i++)
		BOOST_TEST ( a3[i] == A3.vect()[i] );

	/* Ad */

	A3 = SE3Ook::Ad_star(G1,A2);
	a3 = a2.Ad_star(g1);

	for (i=0; i<6; i++)
		BOOST_TEST ( a3[i] == A3.vect()[i] );
}
