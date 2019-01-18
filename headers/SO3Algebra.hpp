#ifndef DEF_LIE_SO3_ALGEBRA
#define DEF_LIE_SO3_ALGEBRA

#include <string>

#include "SO3Group.hpp"

#include "../libs/eigen/Eigen/Dense"
#include "../libs/eigen/Eigen/Geometry"

namespace Lie {
namespace SO3 {

/**
 * Class Lie algebra \f$\mathfrak{so}(3)\f$ implementation.
 * \tparam T Floating point type used for internal representation of coefficients.
 */
template <typename T>
class Algebra {
protected:
	Eigen::Matrix<T,3,1> m_v;

public:
	Algebra<T> ( ) : m_v(Eigen::Matrix<T,3,1>::Zero()) { }

	Algebra<T> (const Eigen::Matrix<T,3,1>& v) : m_v(v) { }

	/* Group operations */

	/**
	 * \return the group inverse of `*this`.
	 */
	Algebra<T> inverse ( ) const {
		Algebra<T> a(-m_v);
		return a;
	}

	/**
	 * Inverts the group element `*this`.
	 */
	void inverted ( ) {
		*this = this->inverse();
	}

	void operator+= (const Algebra<T>& g) {
		m_v += g.v();
	}

	/**
	 * \return the element representing the group identity for operation `+`.
	 */
	static Algebra<T> Identity ( ) {
		Algebra<T> res(Eigen::Matrix<T,3,1>::Zero());
		return res;
	}
	
	/**
	 * \return the group addition of `*this` and \p g.
	 */
	Algebra<T> operator+ (const Algebra<T>& g) const {
		Algebra<T> res(*this);
		res += g;
		return res;
	}

	/* Vector space operations */

	void operator*= (const T& s) {
		m_v *= s;
	}
	
	/**
	 * \return the dot product of `*this` by the scalar \p s.
	 */
	Algebra<T> operator* (const T& s) const {
		Algebra<T> res(*this);
		res *= s;
		return res;
	}

	/* Lie algebra operations */

	/**
	 * Implements the non-commutative Lie bracket operation \f$[a,b]\f$ where \f$a\f$ is `*this`
	 * and \f$b\f$ is \p g.
	 * \return the bracket operation between `*this` and \p g.
	 */
	Algebra<T> bracket (const Algebra<T>& g) const {
		return Algebra<T>(m_v.cross(g.v()));
	}

	/* Other operations */

	/**
	 * Transforms the size 3 vector \p vec by the algebra element represented by `*this`.
	 * With a skew matrix representation, it implements the operation \f$\widehat\omega v\f$;
	 * with a vector reprentation it implements \f$\omega\wedge v\f$.
	 * \return the result of the action of `*this` on \p vec.
	 */
	Eigen::Matrix<T,3,1> operator* (const Eigen::Matrix<T,3,1>& vec) const {
		return m_v.cross(vec);
	}

	/* Accessors */

	Eigen::Quaternion<T> v ( ) const {
		return m_v;
	}

	void v (const Eigen::AngleAxis<T>& aa) {
		m_v = aa.angle()*aa.axis();
	}

	void v (const Eigen::Matrix<T,3,1>& vec) {
		m_v = vec;
	}
	
	// Les fonctions de normalisation sont-elles vraiment utiles ?
	void normalize ( ) {
		m_v.normalize();
	}

	Algebra<T> normalized ( ) const {
		Algebra<T> res(*this);
		res.normalize();
		return res;
	}

	T norm ( ) const {
		return m_v.norm();
	}

	T angle ( ) const {
		return m_v.norm();
	}

	Eigen::Matrix<T,3,1> axis ( ) const {
		return m_v.normalized();
	}

	/**
	 * Implements the exponential map \f$\exp:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 * \return the \ref Lie::SO3::Group<T> implementation of the Lie group \f$SO(3)\f$ element
	 * that represents the exponential of `*this`.
	 */
	Group<T> exp ( ) const {
		T a = m_v.norm()/2.0;
		Eigen::Matrix<T,4,1> V;
		V << cos(a) << sin(a)*m_v.normalized();
		return Group<T>(V);
	}

	/**
	 * Implements the Cayley map \f$cay:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 * \return the \ref Lie::SO3::Group<T> implementation of the Lie group \f$SO(3)\f$ element
	 * that represents the image of `*this` by the Cayley map.
	 */
	Group<T> cay ( ) const {
		T	n =		m_v.norm(),
			den =	4.0+n*n;
		Eigen::Matrix<T,4,1> V;
		V << 1.0 - 2.0*n*n/den << m_v.normalized() * 4.0*n/den;
		return Group<T>(V);
	}

	/**
	 * \return the axis-angle representation of `*this`.
	 */
	Eigen::AngleAxis<T> toAngleAxis ( ) const {
		return Eigen::AngleAxis<T>(m_v.norm(),m_v.normalized());
	}

	Eigen::Matrix<T,3,3> toRotationMatrix ( ) const {
		Eigen::Matrix<T,3,3> mat;
		mat <<	T(0),		-m_v[2],	m_v[1],
				m_v[2],		T(0),		-m_v[0],
				-m_v[1],	m_v[0],		T(0);
		return mat;
	}

	Eigen::Matrix<T,3,1> toVector ( ) const {
		return m_v;
	}

	/**
	 * Implements the three generators of the skew matrix representation of \f$\mathfrak{so}(3)\f$.
	 * Those generators are repsectively
	 * \f[
	 *		\widehat J_0=\begin{pmatrix}0&0&0\\0&0&-1\\0&1&0\end{pmatrix},\quad
	 *		\widehat J_1=\begin{pmatrix}0&0&1\\0&0&0\\-1&0&0\end{pmatrix},\quad
	 *		\widehat J_2=\begin{pmatrix}0&-1&0\\1&0&0\\0&0&0\end{pmatrix}
	 *	\f]
	 *	and define the isomorphism between 3 by 3 skew matrices and size 3 vectors
	 *	\f[ \widehat\omega = \omega\cdot\left(\widehat J_0,\widehat J_1,\widehat J_2\right)^T. \f]
	 * \return the `i`-th generator matrix.
	 */
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

} // namespace SO3
} // namespace Lie

template <typename T>
std::ostream operator<< (std::ostream& stream, Lie::SO3::Algebra<T> const& g) {
	stream << g.q();
	return stream;
}

#endif
