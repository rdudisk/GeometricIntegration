#ifndef DEF_COMMON_ABSTRACT_LIEALGEBRA
#define DEF_COMMON_ABSTRACT_LIEALGEBRA

namespace Abstract
{

template <typename T_SCALAR_TYPE>
class LieAlgebra
{
public:
	static const unsigned int DOF;

public:
	~LieAlgebra<T_SCALAR_TYPE> ( )
	{ }

	/* Group operations */

	/**
	 * \return the inverse of `*this` for the '+' group operation.
	 */
	virtual LieAlgebra<T_SCALAR_TYPE>
	inverse ( ) const = 0;

	/**
	 * Inverts the element `*this` for the '+' group operation.
	 */
	void
	inverted ( )
	{ *this = this->inverse(); }

	virtual void
	operator+= (const LieAlgebra<T_SCALAR_TYPE>& g) = 0
	
	/**
	 * \return the group addition of `*this` and \p g.
	 */
	LieAlgebra<T_SCALAR_TYPE>
	operator+ (const LieAlgebra<T_SCALAR_TYPE>& g) const
	{
		LieAlgebra<T_SCALAR_TYPE> res(*this);
		res += g;
		return res;
	}

	/**
	 * Performs `*this` + \p g.inverse().
	 */
	LieAlgebra<T_SCALAR_TYPE>
	operator- (const LieAlgebra<T_SCALAR_TYPE>& g) const
	{
		LieAlgebra<T_SCALAR_TYPE> res(*this);
		res += g.inverse();
		return res;
	}

	/**
	 * \return the element representing the group identity for operation `+`.
	 */
	virtual static Algebra<T_SCALAR_TYPE>
	Zero ( ) = 0;

	/* Vector space operations */

	virtual void
	operator*= (T_SCALAR_TYPE s) = 0;
	
	/**
	 * \return the dot product of `*this` by the scalar \p s.
	 */
	LieAlgebra<T_SCALAR_TYPE>
	operator* (T_SCALAR_TYPE s) const
	{
		LieAlgebra<T_SCALAR_TYPE> res(*this);
		res *= s;
		return res;
	}

	/* Lie algebra operations */

	/**
	 * Implements the non-commutative Lie bracket operation \f$[a,b]\f$ where \f$a\f$ is `*this`
	 * and \f$b\f$ is \p g.
	 * \return the bracket operation between `*this` and \p g.
	 */
	virtual Algebra<T_SCALAR_TYPE>
	bracket (const Algebra<T_SCALAR_TYPE>& g) const = 0;

	/* Other operations */

	virtual T_SCALAR_TYPE
	norm () const = 0;

	/**
	 * Implements the exponential map \f$\exp:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 * \return the \ref Lie::SO3::Group<T> implementation of the Lie group \f$SO(3)\f$ element
	 * that represents the exponential of `*this`.
	 */
	virtual Abstract::LieGroup
	exp ( ) const = 0;

	virtual NOXVector<DOF>
	toNOXVector ( ) const = 0;

	static const unsigned int
	dof ()
	{
		return DOF;
	}
};

} // namespace Abstract

template <typename T_SCALAR_TYPE>
const Abstract::LieAlgebra<T_SCALAR_TYPE>
operator* (T_SCALAR_TYPE s, const Abstract::LieAlgebra<T_SCALAR_TYPE> g)
{
	return g*s;
}

#endif
