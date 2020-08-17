#ifndef DEF_LIE_SE3_GROUP
#define DEF_LIE_SE3_GROUP

#include <string>

#include "lie_SO3_group.hpp"
#include <Eigen/Dense>

namespace Lie {
namespace SE3 {

template <typename T>
class Group {
protected:
	SO3::Group<T> R_;
	Eigen::Matrix<T,3,1> r_;

public:
	Group<T> ( )
		: R_(SO3::Group<T>::Identity()), r_(Eigen::Matrix<T,3,1>::Zero())
	{ }

	Group<T> (Eigen::Matrix<T,3,3> R, Eigen::Matrix<T,3,1> r)
		: R_(SO3::Group<T>(R)), r_(r)
	{ }

	Group<T> (SO3::Group<T> g, Eigen::Matrix<T,3,1> r)
		: R_(g), r_(r)
	{ }

	Group<T> invert ( ) const {
		Group<T> g(R_.invert(),T(-1)*R_.invert().rotationMatrix()*r_);
		return g;
	}

	void invertInPlace ( ) {
		*this=this->invert;
	}

	void operator*= (Group<T> const& g) {
		r_ += R_.rotationMatrix()*g.translationVect();
		R_ *= g.getSO3();
	}

	Group<T> operator* (Group<T> const& g) const {
		Group<T> res(*this);
		res *= g;
		return res;
	}

	Group<T> operator* (T const& scal) const {
		Group<T> res(*this);
	}

	Eigen::Matrix<T,3,1> operator* (Eigen::Matrix<T,3,1> const& v) const {
		Eigen::Matrix<T,4,1> V;
		V.head(3) = v;
		V(3) = 1;
		return matrixUp()*V;
	}

	SO3::Group<T> getSO3 ( ) const { return R_; }

	void setRotation (T const a, T const b, T const c,\
						std::string const& euler) {
		R_.setRotation(a,b,c,euler);
	}

	void setTranslation (Eigen::Matrix<T,3,1> const& r) {
		r_=r;
	}

	Eigen::Matrix<T,3,3> rotationMatrix ( ) const { return R_.rotationMatrix(); }

	Eigen::Matrix<T,3,1> rotationVect ( ) const { return R_.rotationVect(std::string("zxz")); }
		// Using z-x-z euler convention

	Eigen::Matrix<T,3,1> translationVect ( ) const { return r_; }

	Eigen::Matrix<T,4,4> matrix ( ) const {
		Eigen::Matrix<T,4,4> M;
		M.block(0,0,3,3) = R_.rotationMatrix();
		M.block(0,3,3,1) = r_;
		M.block(3,0,1,4) << 0,0,0,1;
		return M;
	}

	Eigen::Matrix<T,3,4> matrixUp ( ) const {
		Eigen::Matrix<T,4,4> M = (*this).matrix();
		return M.block(0,0,3,4);
	}

	Eigen::Matrix<T,6,1> vect ( ) const {
		Eigen::Matrix<T,6,1> V;
		V.head(3) = rotationVect();
		V.tail(3) = translationVect();
		return V;
	}

	static Group<T> Identity ( ) {
		Group<T> g(SO3::Group<T>::Identity(),Eigen::Matrix<T,3,1>::Zero());
		return g;
	}
};

} // namespace SE3
} // namespace Lie

template <typename T>
std::ostream operator<< (std::ostream& stream, Lie::SE3::Group<T> const& G) {
	stream << G.matrix();
	return stream;
}

#endif
