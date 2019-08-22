#define BOOST_TEST_MODULE "Geomi Common module test"
#include <boost/test/included/unit_test.hpp>

#include <vector>

#include "Geomi/Common"

BOOST_AUTO_TEST_CASE ( common )
{
	int i,j;

	/* Lagrange polynomials tests */

	LagrangePolynomials lp;

	std::vector<double> dates;
	dates.push_back(0.0);
	dates.push_back(0.5);
	dates.push_back(1.0);

	std::vector<std::vector<double>> poly = lp.polynomials(0,dates);

	BOOST_TEST ( poly.size()==dates.size() );
	BOOST_TEST ( poly[0].size()==1 );
	BOOST_TEST ( poly[0][0]==1.0 );
	BOOST_TEST ( poly[2][0]==1.0 );

	poly = lp.polynomials(2,dates);

	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			if (i==j)
				BOOST_TEST ( poly[i][j]==1.0 );
			else
				BOOST_TEST ( poly[i][j]==0.0 );
		}
	}
}
