#ifndef DEF_VEC
#define DEF_VEC

#include "../libs/eigen/Eigen/Dense"

template <typename T,int N>
class Vec : public Eigen::Matrix<T,N,1> {
public:
	void osef() {
		int i = 0;
		i++;
	}
};

#endif
