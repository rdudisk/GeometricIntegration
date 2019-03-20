#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include "test.hpp"

using namespace std;

int main (int argc, char* argv[]) {
	Params params;
	params.m = 1.0;
	DLS dls;

	dls.baselinstep(0.0,0.1,20);
	dls.params(params);

	Q pos0;
	TQ vel0;
	pos0 << -5.0,10.0;
	vel0 << 0.8, 0.0;

	for (int i=0; i<dls.size(); i++) {
		dls.pos(i,pos0);
		dls.vel(i,vel0);
	}

	pos0 << -4.2,10.0;
	dls.pos(1,pos0);

	//dls.evolve(step);
	
	VI vi;
	vi.m_syst = &dls;
	vi.midpoint_fdf();

	dls.write2csv("test.csv");

	return 0;
}
