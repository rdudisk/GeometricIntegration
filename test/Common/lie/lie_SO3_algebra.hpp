#ifndef DEF_LIE_so3
#define DEF_LIE_so3

#include <Eigen/Dense>

namespace Lie {
namespace SO3 {

template <typename T>
class Algebra {
protected:
	Eigen::Matrix<T,3,1> w_;

public:
	Algebra<T> ( ) : w_(Eigen::Matrix<T,3,1>::Zero()) { }

	Algebra<T> (Eigen::Matrix<T,3,1> w) : w_(w) { }

	Algebra<T> (Eigen::Matrix<T,3,3> hatw) {
		w_ << hatw(2,1), hatw(0,2), hatw(1,0);
	}

	Algebra<T> invert ( ) const {
		Eigen::Matrix<T,3,1> W(T(-1)*w_);
		return Algebra<T>(W);
	}

	void invertInPlace ( ) { *this = this->invert(); }

	void operator+= (Algebra<T> const& g) { w_ += g.vect(); }

	Algebra<T> operator+ (Algebra<T> const& g) const {
		Algebra<T> res(*this);
		res += g;
		return res;
	}

	void operator-= (Algebra<T> const& g) { w_ -= g.vect(); }

	Algebra<T> operator- (Algebra<T> const& g) const {
		Algebra<T> res(*this);
		res -= g;
		return res;
	}

	void operator*= (T const& d) { w_ *= d; }

	Algebra<T> operator* (T const& d) const {
		Algebra<T> res(*this);
		res *= d;
		return res;
	}

	void operator/= (T const& d) { w_ /= d; }

	Algebra<T> operator/ (T const& d) const {
		Algebra<T> res(*this);
		res /= d;
		return res;
	}

	Eigen::Matrix<T,3,1> vect ( ) const { return w_; }

	Eigen::Matrix<T,3,3> matrix ( ) const {
		Eigen::Matrix<T,3,3> M;
		M << T(0),-w_(2),w_(1),w_(2),T(0),-w_(0),-w_(1),w_(0),T(0);
		return M;
	}

	Eigen::Matrix<T,3,3> rotationMatrix ( ) const {
		return matrix();
	}

	static Algebra<T> Identity ( ) {
		Eigen::Matrix<T,3,1> V(Eigen::Matrix<T,3,1>::Zero());
		Algebra<T> res(V);
		return res;
	}

	static Eigen::Matrix<T,3,3> GeneratorMatrix (int const i) {
		Eigen::Matrix<T,3,3> res(Eigen::Matrix<T,3,3>::Zero());
		if (i == 0)			{ res(1,2) = T(-1); res(2,1) = T(1);  }
		else if (i == 1)	{ res(0,2) = T(1);  res(2,0) = T(-1); }
		else				{ res(0,1) = T(-1); res(1,0) = T(1);  }
		return res;
	}

	static Eigen::Matrix<T,3,1> GeneratorVect (int const i) {
		Eigen::Matrix<T,3,1> res(Eigen::Matrix<T,3,1>::Zero());
		res(i) = T(1);
		return res;
	}
};

} // namespace Algebra
} // namespace Lie

#endif
