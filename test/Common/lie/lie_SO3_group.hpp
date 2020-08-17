#ifndef DEF_LIE_SO3GROUP
#define DEF_LIE_SO3GROUP

#include <cmath>

#include <Eigen/Dense>

namespace Lie {
namespace SO3 {

template <typename T>
class Group {
protected:
	Eigen::Matrix<T,3,3> R_;

public:
	Group<T> ( ) : R_(Eigen::Matrix<T,3,3>::Identity())
	{ }

	Group<T> (Eigen::Matrix<T,3,3> R) : R_(R)	{
		this->normalise();
	}

	Group<T> invert ( ) const {
		Group<T> g(R_.transpose());
		return g;
	}

	void invertInPlace ( ) {
		*this = this->invert();
	}

	void operator*= (Group<T> const& g) {
		R_ *= g.matrix();
		(*this).normalise();
	}

	Group<T> operator* (Group<T> const& g) const {
		Group<T> res(*this);
		res *= g;
		return res;
	}

	Eigen::Matrix<T,3,1> operator* (Eigen::Matrix<T,3,1> v) const {
		return R_*v;
	}

	Eigen::Matrix<T,3,3> rotationMatrix ( ) const {
		return R_;
	}

	Eigen::Matrix<T,3,3> matrix() const {
		return R_;
	}

	void setRotation (T const a, T const b, T const c, \
							std::string const& euler) {
		R_ = Eigen::Matrix<T,3,3>::Identity();
		T v[3] = {a,b,c};
		for (int i=0; i<3; i++) R_ *= buildEuler(v[i], euler[i]);
		normalise();
	}

	Eigen::Matrix<T,3,3> buildEuler (T const a, char const e) {
		Eigen::Matrix<T,3,3> M(Eigen::Matrix<T,3,3>::Identity());
		T c=cos(a);
		T s=sin(a);
		switch (e) {
		case 'x':
			M << 1,0,0,0,c,-s,0,s,c;
			break;
		case 'y':
			M << c,0,s,0,1,0,-s,0,c;
			break;
		case 'z':
			M << c,-s,0,s,c,0,0,0,1;
			break;
		default:
			break;
		}
		return M;
	}

	Eigen::Matrix<T,3,1> rotationVect (std::string const& euler) const {
		// Warning : only 'z-x-z' convention supported
		// See Wikipedia for the criptic formula
		T Y3 = R_(2,1), Z2 = R_(1,1), Z3 = R_(2,2);
		Eigen::Matrix<T,3,1> res;
		res << acos(Z2/sqrt(1-Z3*Z3)), acos(Z3), acos(Y3/sqrt(1-Z3*Z3));
		return res;
	}

	void normalise ( ) {
		R_ /= cbrt(R_.determinant());
	}

	static Group<T> Identity ( ) {
		Group<T> g(Eigen::Matrix<T,3,3>::Identity());
		return g;
	}
};

} // namespace SO3
} // namespace Lie

#endif
