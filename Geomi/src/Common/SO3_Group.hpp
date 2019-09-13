#ifndef DEF_COMMON_SO3_GROUP
#define DEF_COMMON_SO3_GROUP

#include <string>

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace SO3
{

/**
 * Class for Lie group \f$SO(3)\f$ implementation.
 * \tparam T_SCALAR_TYPE Floating point type used for internal representation of coefficients.
 */
template <typename T_SCALAR_TYPE>
class Group : LieGroupBase<Group<T_SCALAR_TYPE>,3>
{
protected:
	Eigen::Quaternion<T_SCALAR_TYPE> m_q;

//public:
	//static const unsigned int DOF = 3;

public:
	Group<T_SCALAR_TYPE> ( )
	: m_q(Eigen::Quaternion<T_SCALAR_TYPE>::Identity())
	{ }

	Group<T_SCALAR_TYPE> (const Group<T_SCALAR_TYPE>& g)
	: m_q(g.m_q)
	{ }

	Group<T_SCALAR_TYPE> (const Eigen::Quaternion<T_SCALAR_TYPE>& q_)
	{
		q(q_);
	}

	Group<T_SCALAR_TYPE>
		(	const T_SCALAR_TYPE& w,
			const T_SCALAR_TYPE& x,
			const T_SCALAR_TYPE& y,
			const T_SCALAR_TYPE& z)
	{
		q(w,x,y,z);
	}

	Group<T_SCALAR_TYPE> (const T_SCALAR_TYPE* data)
	{
		q(data);
	}

	template<class Derived>
	Group<T_SCALAR_TYPE> (const Eigen::QuaternionBase<Derived>& other)
	{
		q(other);
	}

	Group<T_SCALAR_TYPE> (const Eigen::AngleAxis<T_SCALAR_TYPE>& aa)
	{
		q(aa);
	}

	template<typename Derived>
	Group<T_SCALAR_TYPE> (const Eigen::MatrixBase<Derived>& other)
	{
		q(other);
	}

	template<class OtherScalar, int OtherOptions>
	Group<T_SCALAR_TYPE> (const Eigen::Quaternion<OtherScalar,OtherOptions>& other)
	{
		q(other);
	}

	/* Accessors and setters */

	/**
	 * \return the quaternion representation of the rotation.
	 */
	Eigen::Quaternion<T_SCALAR_TYPE>
	q ( ) const
	{
		return m_q;
	}

	void
	q (const Eigen::Quaternion<T_SCALAR_TYPE>& q_)
	{
		m_q = q_;
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	void
	q (	const T_SCALAR_TYPE& w,
		const T_SCALAR_TYPE& x,
		const T_SCALAR_TYPE& y,
		const T_SCALAR_TYPE& z)
	{
		m_q = Eigen::Quaternion<T_SCALAR_TYPE>(w,x,y,z);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	void
	q (const T_SCALAR_TYPE* data)
	{
		m_q = Eigen::Quaternion<T_SCALAR_TYPE>(data);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	template<class Derived>
	void
	q (const Eigen::QuaternionBase<Derived>& other)
	{
		m_q = Eigen::Quaternion<T_SCALAR_TYPE>(other);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	void
	q (const Eigen::AngleAxis<T_SCALAR_TYPE>& aa)
	{
		m_q = Eigen::Quaternion<T_SCALAR_TYPE>(aa);
		m_q.normalize();
	}
	
	// see Eigen::Quaternion constructors
	template<typename Derived>
	void
	q (const Eigen::MatrixBase<Derived>& other)
	{
		m_q = Eigen::Quaternion<T_SCALAR_TYPE>(other);
		m_q.normalize();
	}

	// see Eigen::Quaternion constructors
	template<class OtherScalar, int OtherOptions>
	void
	q (const Eigen::Quaternion<OtherScalar,OtherOptions>& other)
	{
		m_q = Eigen::Quaternion<T_SCALAR_TYPE>(other);
		m_q.normalize();
	}

	/* Group operations */

	/**
	 * @return the group inverse of `*this`.
	 */
	Group<T_SCALAR_TYPE>
	inverse( ) const
	{
		Group<T_SCALAR_TYPE> g(m_q.inverse());
		return g;
	}

	/**
	 * Inverts the group element `*this`
	 */
	/*void
	inverted ( )
	{
		*this=this->invert;
	}*/

	/**
	 * \return the element representing the group identity for operation `*`.
	 */
	static Group<T_SCALAR_TYPE>
	Identity ( )
	{
		Group<T_SCALAR_TYPE> g(Eigen::Quaternion<T_SCALAR_TYPE>::Identity());
		return g;
	}

	/* Group operation is '*' and not '+' since we are used to the matrix representation,
	 * in which case the mutliplication is the group operation */
	void
	operator*= (Group<T_SCALAR_TYPE> const& g)
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
	Group<T_SCALAR_TYPE>
	operator* (Group<T_SCALAR_TYPE> const& g) const
	{
		Group<T_SCALAR_TYPE> res(*this);
		res *= g;
		return res;
	}

	/**
	 * \return the transformation of the \f$\mathbb R^3\f$ vector \p v by the rotation
	 * represented by `*this`.
	 */
	Eigen::Matrix<T_SCALAR_TYPE,3,1>
	transformVector (const Eigen::Matrix<T_SCALAR_TYPE,3,1>& v) const
	{
		return this->toRotationMatrix()*v;
	}

	/**
	 * Implements the product of the matrix representation of the rotation `*this`
	 * and the voctor \p v.
	 * \return the transformation of the \f$\mathbb R^3\f$ vector \p v by the rotation
	 * represented by `*this`.
	 */
	Eigen::Matrix<T_SCALAR_TYPE,3,1>
	operator* (Eigen::Matrix<T_SCALAR_TYPE,3,1> const& v) const
	{
		return m_q._transformVector(v);
	}

	/**
	 * \return the 3 by 3 matrix representation of the rotation.
	 */
	Eigen::Matrix<T_SCALAR_TYPE,3,3>
	toRotationMatrix ( ) const
	{
		return m_q.toRotationMatrix();
	}

	/**
	 * \return the axis-angle representation of the rotation.
	 */
	Eigen::Matrix<T_SCALAR_TYPE,3,3>
	toAxisAngle ( ) const
	{
		return Eigen::AngleAxis<T_SCALAR_TYPE>(m_q);
	}

	/**
	 * \return the vector representation of the rotation, that is the vector \f$\vec v\f$
	 * such that the rotation of any given vector \f$\vec u\f$ is the result of \f$vec v\wedge\vec u\f$.
	 */
	Eigen::Matrix<T_SCALAR_TYPE,3,1>
	toVector ( ) const
	{
		Eigen::AngleAxis<T_SCALAR_TYPE> aa(m_q);
		return aa.angle()*aa.axis();
	}

	bool
	isApprox (Group<T_SCALAR_TYPE> const& g) const
	{
		return m_q.isApprox(g.q());
	}

	/*static const unsigned int
	dof ()
	{
		return DOF;
	}*/
};

} // namespace SO3

template <typename T_SCALAR_TYPE>
std::ostream
operator<< (std::ostream& stream, SO3::Group<T_SCALAR_TYPE> const& G)
{
	stream << G.toRotationMatrix();
	return stream;
}

#endif
