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
	dls.baselinstep(0.0,0.1,400);
	Q pos0;
	TQ vel0;
	pos0 << -5.0,100.0;
	vel0 << 0.8, 0.0;
	//dls.pos(0,Vec2((float)-5.0,(float)100.0));
	//dls.vel(0,Vec2(0.8,0.0));
	dls.pos(0,pos0);
	dls.vel(0,vel0);
	dls.params(params);
	//cout << dls << endl;
	dls.evolve(step);
	//cout << dls << endl;
	dls.write2csv("test.csv");
	return 0;
}
