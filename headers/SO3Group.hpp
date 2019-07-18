#ifndef DEF_LIE_SO3_GROUP
#define DEF_LIE_SO3_GROUP

#include <string>

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace Lie
{
namespace SO3
{

/**
 * Class for Lie group \f$SO(3)\f$ implementation.
 * \tparam T Floating point type used for internal representation of coefficients.
 */
template <typename T>
class Group {
protected:
	Eigen::Quaternion<T> m_q;

public:
	Group<T> ( )
	: m_q(Eigen::Quaternion<T>::Identity())
	{ }

	Group<T> (const Group<T>& g)
	: m_q(g.m_q)
	{ }

	Group<T> (const Eigen::Quaternion<T>& q_)
	{
		q(q_);
	}

	Group<T> (const T& w, const T& x, const T& y, const T& z)
	{
		q(w,x,y,z);
	}

	Group<T> (const T* data)
	{
		q(data);
	}

	template<class Derived>
	Group<T> (const Eigen::QuaternionBase<Derived>& other)
	{
		q(other);
	}

	Group<T> (const Eigen::AngleAxis<T>& aa)
	{
		q(aa);
	}

	template<typename Derived>
	Group<T> (const Eigen::MatrixBase<Derived>& other)
	{
		q(other);
	}

	template<class OtherScalar, int OtherOptions>
	Group<T> (const Eigen::Quaternion<OtherScalar,OtherOptions>& other)
	{
		q(other);
	}

	/* Accessors and setters */

	/**
	 * \return the quaternion representation of the rotation.
	 */
	Eigen::Quaternion<T>
	q ( ) const
	{
		return m_q;
	}

	void
	q (const Eigen::Quaternion<T>& q_)
	{
		m_q = q_;
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	void
	q (const T& w, const T& x, const T& y, const T& z)
	{
		m_q = Eigen::Quaternion<T>(w,x,y,z);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	void
	q (const T* data)
	{
		m_q = Eigen::Quaternion<T>(data);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	template<class Derived>
	void
	q (const Eigen::QuaternionBase<Derived>& other)
	{
		m_q = Eigen::Quaternion<T>(other);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	void
	q (const Eigen::AngleAxis<T>& aa)
	{
		m_q = Eigen::Quaternion<T>(aa);
		m_q.normalize();
	}
	
	// see Eigen::Quaternion constructors
	template<typename Derived>
	void
	q (const Eigen::MatrixBase<Derived>& other)
	{
		m_q = Eigen::Quaternion<T>(other);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	template<class OtherScalar, int OtherOptions>
	void
	q (const Eigen::Quaternion<OtherScalar,OtherOptions>& other)
	{
		m_q = Eigen::Quaternion<T>(other);
		m_q.normalize();
	}

	static size_t
	dof ( )
	{
		return 3;
	}

	/* Group operations */

	/**
	 * @return the group inverse of `*this`.
	 */
	Group<T>
	inverse( ) const
	{
		Group<T> g(m_q.inverse());
		return g;
	}

	/**
	 * Inverts the group element `*this`
	 */
	void
	inverted ( )
	{
		*this=this->invert;
	}

	/**
	 * \return the element representing the group identity for operation `*`.
	 */
	static Group<T>
	Identity ( )
	{
		Group<T> g(Eigen::Quaternion<T>::Identity());
		return g;
	}

	/* Group operation is '*' and not '+' since we are used to the matrix representation,
	 * in which case the mutliplication is the group operation */
	void
	operator*= (Group<T> const& g)
	{
		m_q *= g.q();
		m_q.normalize();
	}

	/**
	 * Implements the group operation on \f$SO(3)\f$.
	 * The `*` operator has been chosen instead of the `+` operator since the usual
	 * matrix and quaternion representations of rotations are groups defined with a
	 * natural mutliplication operation.
	 * \return the group addition of `*this` and \p g.
	 */
	Group<T>
	operator* (Group<T> const& g) const
	{
		Group<T> res(*this);
		res *= g;
		return res;
	}

	/**
	 * \return the transformation of the \f$\mathbb R^3\f$ vector \p v by the rotation
	 * represented by `*this`.
	 */
	Eigen::Matrix<T,3,1>
	transformVector (const Eigen::Matrix<T,3,1>& v) const
	{
		return this->toRotationMatrix()*v;
	}

	/**
	 * Implements the product of the matrix representation of the rotation `*this`
	 * and the voctor \p v.
	 * \return the transformation of the \f$\mathbb R^3\f$ vector \p v by the rotation
	 * represented by `*this`.
	 */
	Eigen::Matrix<T,3,1>
	operator* (Eigen::Matrix<T,3,1> const& v) const
	{
		return m_q.transformVector(v);
	}

	/**
	 * \return the 3 by 3 matrix representation of the rotation.
	 */
	Eigen::Matrix<T,3,3>
	toRotationMatrix ( ) const
	{
		return m_q.toRotationMatrix();
	}

	/**
	 * \return the axis-angle representation of the rotation.
	 */
	Eigen::Matrix<T,3,3>
	toAxisAngle ( ) const
	{
		return Eigen::AngleAxis<T>(m_q);
	}

	/**
	 * \return the vector representation of the rotation, that is the vector \f$\vec v\f$
	 * such that the rotation of any given vector \f$\vec u\f$ is the result of \f$vec v\wedge\vec u\f$.
	 */
	Eigen::Matrix<T,3,1>
	toVector ( ) const
	{
		Eigen::AngleAxis<T> aa(m_q);
		return aa.angle()*aa.axis();
	}

	bool
	isApprox (Group<T> const& g) const
	{
		return m_q.isApprox(g.q());
	}
};

} // namespace SO3
} // namespace Lie

template <typename T>
std::ostream
operator<< (std::ostream& stream, Lie::SO3::Group<T> const& G)
{
	stream << G.toRotationMatrix();
	return stream;
}

#endif
