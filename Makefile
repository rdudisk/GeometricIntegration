all: reissner

reissner: main.cpp headers/*
	#g++ -o main main.cpp -lgsl -lgslcblas -lm -lpthread -std=gnu++11 
	g++ -o main main.cpp libs/fastgl/fastgl.cpp -lgsl -lgslcblas -lm  -std=c++11

test: test-so3.cpp headers/*
	g++ -o test test-so3.cpp -lgsl -lgslcblas -lm -lpthread -std=c++11

doxygen:
	doxygen Doxygen
