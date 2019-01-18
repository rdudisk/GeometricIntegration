#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include "headers/RigidBody.hpp"

using namespace std;

int main (int argc, char* argv[]) {
	Params<T> params;
	params.I <<	1, 0, 0,
				0, 2, 0,
				0, 0, 4;
	DLS dls;
	dls.baselinstep(0.0,0.1,400);
	dls.params(params);
	dls.evolve_exact();
	dls.write2csv("test-so3.csv");
	return 0;
}
