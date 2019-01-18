#ifndef DEF_LIE_SE3_GROUP
#define DEF_LIE_SE3_GROUP

#include <string>

#include "../libs/eigen/Eigen/Dense"
#include "../libs/eigen/Eigen/Geometry"

namespace Lie {
namespace SE3 {

template <typename T>
class Group {
protected:
	Eigen::Quaternion<T> m_q;
	Eigen::Matrix<T,3,1> m_r;

public:
	Group<T> ( )
		: m_q(Eigen::Quaternion<T>::Identity()), m_r(Eigen::Matrix<T,3,1>::Zero())
	{ }

	Group<T> (const Group<T>& g)
		: m_q(g.m_q), m_r(g.m_r)
	{ }

	Group<T> (Eigen::Quaternion<T> q, Eigen::Matrix<T,3,1> r)
		: m_q(q), m_r(r)
	{ }

	~Group<T> ( ) {
		delete m_q;
		delete m_r;
	}

	Eigen::Quaternion<T> q ( ) const {
		return m_q;
	}

	Eigen::Matrix<T,3,1> r ( ) const {
		return m_r;
	}

	void q (const Eigen::Quaternion<T>& q_) {
		m_q = q_;
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	void q (const T& w, const T& x, const T& y, const T& z) {
		q = Eigen::Quaternion<T>(w,x,y,z);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	void q (const T* data) {
		q = Eigen::Quaternion<T>(data);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	template<class Derived>
	void q (const Eigen::QuaternionBase<Derived>& other) {
		q = Eigen::Quaternion<T>(other);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	void q (const Eigen::AngleAxis<T>& aa) {
		q = Eigen::Quaternion<T>(aa);
		m_q.normalize();
	}
	
	// see Eigen::Quaternion constructors
	template<class Derived>
	void q (const Eigen::MatrixBase<Derived>& other) {
		q = Eigen::Quaternion<T>(other);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	template<class OtherScalar, int OtherOptions>
	void q (const Eigen::Quaternion<OtherScalar,OtherOptions>& other) {
		q = Eigen::Quaternion<T>(other);
		m_q.normalize();
	}

	void r (const Eigen::Matrix<T,3,1>& r_) {
		m_r = r_;
	}

	Group<T> inverse( ) const {
		Group<T> g(m_q.inverse(),m_q.inverse().transformVector(m_r));
		return g;
	}

	void inverted ( ) {
		*this=this->invert;
	}

	void operator*= (Group<T> const& g) {
		m_r += m_q.transformVector(g.r());
		m_q *= g.q();
		m_q.normalize();
	}

	Group<T> operator* (Group<T> const& g) const {
		Group<T> res(*this);
		res *= g;
		return res;
	}

	/*
	Group<T> operator* (T const& scal) const {
		Group<T> res(*this);
		// ???
		return res;
	}
	*/

	Eigen::Matrix<T,3,1> transformVector (const Eigen::Matrix<T,3,1>& v) const {
		return (*this)*v;
	}

	Eigen::Matrix<T,3,1> operator* (Eigen::Matrix<T,3,1> const& v) const {
		return m_q.transformVector(v)+m_r;
	}

	/*
	void setRotation (T const a, T const b, T const c,\
						std::string const& euler) {
		m_q.setRotation(a,b,c,euler);
	}
	*/

	Eigen::Matrix<T,3,3> toRotationMatrix ( ) const {
		return m_q.toRotationMatrix();
	}

	/*
	Eigen::Matrix<T,3,1> toRotationVector ( ) const {
		// Using z-x-z euler convention
		return m_q.rotationVect(std::string("zxz"));
	}
	*/

	Eigen::Matrix<T,3,1> toTranslationVector ( ) const {
		return m_r;
	}

	Eigen::Matrix<T,4,4> matrix ( ) const {
		Eigen::Matrix<T,4,4> M;
		M.block(0,0,3,3) = m_q.toRotationMatrix();
		M.block(0,3,3,1) = m_r;
		M.block(3,0,1,4) << 0,0,0,1;
		return M;
	}

	/*
	Eigen::Matrix<T,3,4> matrixUp ( ) const {
		Eigen::Matrix<T,4,4> M = (*this).matrix();
		return M.block(0,0,3,4);
	}
	*/

	/*
	Eigen::Matrix<T,6,1> vector ( ) const {
		Eigen::Matrix<T,6,1> V;
		V.head(3) = toRotationVector();
		V.tail(3) = m_r;
		return V;
	}
	*/

	static Group<T> Identity ( ) {
		Group<T> g(Eigen::Quaternion<T>::Identity(),Eigen::Matrix<T,3,1>::Zero());
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
