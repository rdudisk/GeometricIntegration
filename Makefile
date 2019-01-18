all: main

reissner: main.cpp headers/*
	#g++ -o main main.cpp -lgsl -lgslcblas -lm -lpthread -std=gnu++11 
	g++ -o main main.cpp problem1.cpp -lgsl -lgslcblas -lm -lpthread -std=c++11 

test: test-so3.cpp headers/*
	g++ -o test test-so3.cpp -lgsl -lgslcblas -lm -lpthread -std=c++11
