#ifndef DEF_LIE_se3
#define DEF_LIE_se3

#include <Eigen/Dense>
#include "lie_SO3_algebra.hpp"

using namespace Eigen;

namespace Lie {
namespace SE3 {

template <typename T>
class Algebra {
protected:
	SO3::Algebra<T> w_;
	Eigen::Matrix<T,3,1> gam_;

public:
	Algebra<T> ( )
		: w_(SO3::Algebra<T>::Identity()), gam_(Eigen::Matrix<T,3,1>::Zero())
	{ }

	Algebra<T> (Eigen::Matrix<T,3,1> w, Eigen::Matrix<T,3,1> gam)
		: w_(w), gam_(gam)
	{ }

	Algebra<T> (SO3::Algebra<T> g, Eigen::Matrix<T,3,1> gam)
		: w_(g), gam_(gam)
	{ }

	Algebra<T> (Eigen::Matrix<T,6,1> V)
		: w_(Eigen::Matrix<T,3,1>(V.head(3))) ,gam_(V.tail(3))
	{ }

	Algebra<T> invert ( ) const {
		return Algebra<T>(w_.invert(),T(-1)*gam_);
	}

	Algebra<T> invertInPlace ( ) const {
		*this = this->invert();
	}

	void operator+= (Algebra<T> const& g) {
		w_ += g.getSO3();
		gam_ += g.translationVect();
	}

	Algebra<T> operator+ (Algebra<T> const& g) const {
		Algebra<T> res(*this);
		res += g;
		return res;
	}

	void operator-= (Algebra<T> const& g) {
		w_ -= g.getSO3();
		gam_ -= g.translationVect();
	}

	Algebra<T> operator- (Algebra<T> const& g) const {
		Algebra<T> res(*this);
		res -= g;
		return res;
	}

	void operator*= (T const& d) {
		w_ *= d;
		gam_ *= d;
	}

	Algebra<T> operator* (T const& d) const {
		Algebra<T> res(*this);
		res *= d;
		return res;
	}

	void operator/= (T const& d) {
		w_ /= d;
		gam_ /= d;
	}

	Algebra<T> operator/ (T const& d) const {
		Algebra<T> res(*this);
		res /= d;
		return res;
	}
	
	Eigen::Matrix<T,6,1> vect ( ) const {
		Eigen::Matrix<T,6,1> V;
		V.head(3) = w_.vect();
		V.tail(3) = gam_;
		return V;
	}

	Eigen::Matrix<T,4,4> matrix ( ) const {
		Eigen::Matrix<T,4,4> M(Eigen::Matrix<T,4,4>::Zero());
		M.block(0,0,3,3) = w_.matrix();
		M.block(0,3,3,1) = gam_;
		return M;
	}

	SO3::Algebra<T> getSO3 ( ) const {
		return w_;
	}

	Eigen::Matrix<T,3,1> rotationVect ( ) const {
		return w_.vect();
	}

	Eigen::Matrix<T,3,3> rotationMatrix ( ) const {
		return w_.matrix();
	}

	Eigen::Matrix<T,3,1> translationVect ( ) const {
		return gam_;
	}

	static Algebra<T> Identity ( ) {
		Algebra<T> g(SO3::Algebra<T>::Identity(),Eigen::Matrix<T,3,1>::Zero());
		return g;
	}

	static Eigen::Matrix<T,4,4> GeneratorMatrix (int const i) {
		Eigen::Matrix<T,4,4> res(Eigen::Matrix<T,4,4>::Zero());
		if (i<3)	res.block(0,0,3,3) = SO3::Algebra<T>::GeneratorMatrix(i);
		else		res(i-3,3) = T(1);
		return res;
	}

	static Eigen::Matrix<T,6,1> GeneratorVect (int const i) {
		Eigen::Matrix<T,6,1> res(Eigen::Matrix<T,6,1>::Zero());
		if (i<3)	res.head(3) = SO3::Algebra<T>::GeneratorVect(i);
		else		res(i,0) = T(1);
		return res;
	}
};

} // namespace SE3
} // namespace Lie

template <typename T>
std::ostream operator<< (std::ostream& stream, Lie::SE3::Algebra<T> const& g) {
	stream << g.matrix();
	return stream;
}

#endif
