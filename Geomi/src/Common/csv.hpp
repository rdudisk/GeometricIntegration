#ifndef DEF_COMMON_CSV
#define DEF_COMMON_CSV

#include <Eigen/Dense>

//#include "../../Space/SO3Group.hpp"
//#include "../../Space/SO3Algebra.hpp"
//#include "Common_NOXVector.hpp"

template <typename T>
std::string
csvString (const T f, const std::string sep);

template <>
std::string
csvString (const float f, const std::string sep)
{
	std::ostringstream ss;
	ss << f;
	std::string str(ss.str());
	return str;
}

template <>
std::string
csvString (const double f, const std::string sep)
{
	std::ostringstream ss;
	ss << f;
	std::string str(ss.str());
	return str;
}

template <int N>
std::string
csvString (const Eigen::Matrix<float,N,1> vec, const std::string sep)
{
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
std::string
csvString (const Eigen::Matrix<float,2,1> vec, const std::string sep)
{
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
std::string
csvString (const SO3::Group<float> q, const std::string sep)
{
	std::ostringstream ss;
	Eigen::Matrix<float,3,3> m = q.rotationMatrix();
	for (int i=0; i<9; i++) {
		ss << m(i/3,i%3);
		if (i!=8)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string
csvString (const SO3::Algebra<float> q, const std::string sep)
{
	std::ostringstream ss;
	Eigen::Matrix<float,3,1> v = q.vector();
	for (int i=0; i<3; i++) {
		ss << v[i];
		if (i!=2)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string
csvString (const SO3::Group<double> q, const std::string sep)
{
	std::ostringstream ss;
	Eigen::Matrix<double,3,3> m = q.rotationMatrix();
	for (int i=0; i<9; i++) {
		ss << m(i/3,i%3);
		if (i!=8)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string
csvString (const SO3::Algebra<double> q, const std::string sep)
{
	std::ostringstream ss;
	Eigen::Matrix<double,3,1> v = q.vector();
	for (int i=0; i<3; i++) {
		ss << v[i];
		if (i!=2)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <int N>
std::string
csvString (const NOXVector<N> vec, const std::string sep)
{
	std::ostringstream ss;
	for (int i=0,r=vec.size()-1; i<=r; i++) {
		ss << vec[i];
		if (i<r)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string
csvString (const NOXVector<2> vec, const std::string sep)
{
	std::ostringstream ss;
	for (int i=0,r=2; i<=r; i++) {
		ss << vec[i];
		if (i<r)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string
csvString (const NOXVector<3> vec, const std::string sep)
{
	std::ostringstream ss;
	for (int i=0,r=2; i<=r; i++) {
		ss << vec[i];
		if (i<r)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

template <>
std::string
csvString (const NOXVector<6> vec, const std::string sep)
{
	std::ostringstream ss;
	for (int i=0,r=5; i<=r; i++) {
		ss << vec[i];
		if (i<r)
			ss << sep;
	}
	std::string str(ss.str());
	return str;
}

#endif
