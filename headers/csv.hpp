#ifndef DEF_CSV
#define DEF_CSV

#include "../libs/eigen/Eigen/Dense"
#include "Vec.hpp"
#include "SO3Group.hpp"
#include "SO3Algebra.hpp"

template <typename T>
std::string csvString(const T f, const std::string sep);

template <>
std::string csvString(const float f, const std::string sep) {
	std::ostringstream ss;
	ss << f;
	std::string str(ss.str());
	return str;
}

template <int N>
std::string csvString(const Eigen::Matrix<float,N,1> vec, const std::string sep) {
	std::ostringstream ss;
	for (int i=0,r=vec.rows()-1; i<=r; i++) {
		ss << vec[i];
		if (i<r)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string csvString(const Eigen::Matrix<float,2,1> vec, const std::string sep) {
	std::ostringstream ss;
	for (int i=0; i<=1; i++) {
		ss << vec[i];
		if (i<1)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string csvString(const Vec<float,2> vec, const std::string sep) {
	std::ostringstream ss;
	for (int i=0; i<=1; i++) {
		ss << vec[i];
		if (i<1)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string csvString(const Lie::SO3::Group<float> q, const std::string sep) {
	std::ostringstream ss;
	Eigen::Matrix<float,3,3> m = q.toRotationMatrix();
	for (int i=0; i<9; i++) {
		ss << m(i/3,i%3);
		if (i!=8)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string csvString(const Lie::SO3::Algebra<float> q, const std::string sep) {
	std::ostringstream ss;
	Eigen::Matrix<float,3,3> m = q.toRotationMatrix();
	for (int i=0; i<9; i++) {
		ss << m(i/3,i%3);
		if (i!=8)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

#endif
