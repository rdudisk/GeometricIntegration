#ifndef DEF_COMMON_SE3_GROUP
#define DEF_COMMON_SE3_GROUP

#include <string>

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace SE3
{

/**
 * Class for Lie group \f$SO(3)\f$ implementation.
 * \tparam T_SCALAR Floating point type used for internal representation of coefficients.
 */
template <typename T_SCALAR>
class Group : public LieGroupBase<Group<T_SCALAR>,6>
{
protected:
	Eigen::Quaternion<T_SCALAR> m_q;
	Eigen::Matrix<T_SCALAR,3,1> m_trans;

//public:
	//static const unsigned int DOF = 3;

public:
	Group ( )
	: m_q(Eigen::Quaternion<T_SCALAR>::Identity()), m_trans(Eigen::Matrix<T_SCALAR,3,1>::Zero())
	{ }

	Group (const Group<T_SCALAR>& g)
	: m_q(g.m_q), m_trans(g.m_trans)
	{ }

	Group (const Eigen::Quaternion<T_SCALAR>& q, const Eigen::Matrix<T_SCALAR,3,1>& trans)
	: m_q(q), m_trans(trans)
	{ }

	Group (const Eigen::Matrix<T_SCALAR,4,4>& m)
	: m_trans(m.block(0,3,3,1))
	{
		Eigen::Matrix<T_SCALAR,3,3> tmp = m.block(0,0,3,3);
		m_q = Eigen::Quaternion<T_SCALAR>(tmp);
	}

	/*
	Group<T_SCALAR> (	const T_SCALAR& w,
							const T_SCALAR& x,
							const T_SCALAR& y,
							const T_SCALAR& z)
	{ q(w,x,y,z); }

	Group (const T_SCALAR* data)
	{ q(data); }

	template<class Derived>
	Group<T_SCALAR> (const Eigen::QuaternionBase<Derived>& other)
	{ q(other); }

	Group<T_SCALAR> (const Eigen::AngleAxis<T_SCALAR>& aa)
	{ q(aa); }

	template<typename Derived>
	Group<T_SCALAR> (const Eigen::MatrixBase<Derived>& other)
	{ q(other); }

	template<class OtherScalar, int OtherOptions>
	Group<T_SCALAR> (const Eigen::Quaternion<OtherScalar,OtherOptions>& other)
	{ q(other); }
	*/

	/* Accessors and setters */

	Eigen::Matrix<T_SCALAR,3,1>
	trans ( ) const
	{ return m_trans; }

	void
	trans (const Eigen::Matrix<T_SCALAR,3,1>& v)
	{ m_trans = v; }
	
	void
	trans (const size_t i, const T_SCALAR s)
	// TODO check index domain
	{ m_trans[i] = s; }

	/**
	 * \return the quaternion representation of the rotation.
	 */
	Eigen::Quaternion<T_SCALAR>
	q ( ) const
	{ return m_q; }

	void
	q (const Eigen::Quaternion<T_SCALAR>& q_)
	{ m_q = q_; m_q.normalize(); }

	// see Eigen::Quaternion constructors
	void
	q (	const T_SCALAR& w,
		const T_SCALAR& x,
		const T_SCALAR& y,
		const T_SCALAR& z)
	{ m_q = Eigen::Quaternion<T_SCALAR>(w,x,y,z); m_q.normalize(); }

	// see Eigen::Quaternion constructors
	void
	q (const T_SCALAR* data)
	{ m_q = Eigen::Quaternion<T_SCALAR>(data); m_q.normalize(); }

	// see Eigen::Quaternion constructors
	template<class Derived>
	void
	q (const Eigen::QuaternionBase<Derived>& other)
	{ m_q = Eigen::Quaternion<T_SCALAR>(other); m_q.normalize(); }

	// see Eigen::Quaternion constructors
	void
	q (const Eigen::AngleAxis<T_SCALAR>& aa)
	{ m_q = Eigen::Quaternion<T_SCALAR>(aa); m_q.normalize(); }
	
	// see Eigen::Quaternion constructors
	template<typename Derived>
	void
	q (const Eigen::MatrixBase<Derived>& other)
	{ m_q = Eigen::Quaternion<T_SCALAR>(other); m_q.normalize(); }

	// see Eigen::Quaternion constructors
	template<class OtherScalar, int OtherOptions>
	void
	q (const Eigen::Quaternion<OtherScalar,OtherOptions>& other)
	{ m_q = Eigen::Quaternion<T_SCALAR>(other); m_q.normalize(); }

	/*
	T_SCALAR const&
	operator[] (size_t index) const
	{ return m_trans[index-3]; }

	T_SCALAR&
	operator[] (size_t index)
	{ return m_trans[index-3]; }
	*/

	/* Group operations */

	/**
	 * @return the group inverse of `*this`.
	 */
	Group<T_SCALAR>
	inverse( ) const
	{ return Group<T_SCALAR>(m_q.inverse(),-m_q.toRotationMatrix().transpose()*m_trans); }

	/**
	 * \return the element representing the group identity for operation `*`.
	 */
	static Group<T_SCALAR>
	Identity ( )
	{ return Group<T_SCALAR>(Eigen::Quaternion<T_SCALAR>::Identity(),Eigen::Matrix<T_SCALAR,3,1>::Zero()); }

	static Group<T_SCALAR>
	Random ( )
	{ return Group<T_SCALAR>(Eigen::Quaternion<T_SCALAR>::UnitRandom(),Eigen::Matrix<T_SCALAR,3,1>::Zero()); }

	/* Group operation is '*' and not '+' since we are used to the matrix representation,
	 * in which case the mutliplication is the group operation */
	void
	operator*= (Group<T_SCALAR> const& g)
	{ m_q *= g.q(); m_q.normalize(); m_trans += g.trans(); }

	/**
	 * Implements the group operation on \f$SO(3)\f$.
	 * The `*` operator has been chosen instead of the `+` operator since the usual
	 * matrix and quaternion representations of rotations are groups defined with a
	 * natural mutliplication operation.
	 * \return the group addition of `*this` and \p g.
	 */
	Group<T_SCALAR>
	operator* (Group<T_SCALAR> const& g) const
	{
		Group<T_SCALAR> res(*this);
		res *= g;
		return res;
	}


	/**
	 * \return the transformation of the \f$\mathbb R^3\f$ vector \p v by `*this`.
	 */
	Eigen::Matrix<T_SCALAR,3,1>
	transformVector (const Eigen::Matrix<T_SCALAR,3,1>& v) const
	{
		return this->rotationMatrix()*v+this->m_trans;
	}

	Eigen::Matrix<T_SCALAR,3,1>
	rotateVector (const Eigen::Matrix<T_SCALAR,3,1>& v) const
	{
		return this->rotationMatrix()*v;
	}

	/**
	 * Implements the product of the matrix representation of the rotation `*this`
	 * and the voctor \p v.
	 * \return the transformation of the \f$\mathbb R^3\f$ vector \p v by the rotation
	 * represented by `*this`.
	 */
	Eigen::Matrix<T_SCALAR,3,1>
	operator* (Eigen::Matrix<T_SCALAR,3,1> const& v) const
	{ return this->transformVector(v); }

	/**
	 * \return the 3 by 3 matrix representation of the rotation.
	 */
	Eigen::Matrix<T_SCALAR,3,3>
	rotationMatrix ( ) const
	{ return m_q.toRotationMatrix(); }

	Eigen::Matrix<T_SCALAR,4,4>
	matrix ( ) const
	{
		Eigen::Matrix<T_SCALAR,4,4> M = Eigen::Matrix<T_SCALAR,4,4>::Zero();
		M.block(0,0,3,3) = this->rotationMatrix();
		M.block(0,3,3,1) = m_trans;
		M(3,3) = T_SCALAR(1);
		return M;
	}

	bool
	isApprox (Group<T_SCALAR> const& g) const
	{ return m_q.isApprox(g.q()) & m_trans.isApprox(g.trans()); }

};

} // namespace SO3

template <typename T_SCALAR>
std::ostream
operator<< (std::ostream& stream, SE3::Group<T_SCALAR> const& G)
{
	stream << G.matrix();
	return stream;
}

#endif
