#ifndef DEF_COMMON_LIEALGEBRABASE
#define DEF_COMMON_LIEALGEBRABASE

template <typename T_DERIVED, typename T_GROUP, unsigned int T_DOF, typename T_SCALAR_TYPE>
class LieAlgebraBase : CRTP<T_DERIVED>
{
public:
	static const unsigned int DOF = T_DOF;

public:
	/* Group operations */

	/**
	 * \return the inverse of `*this` for the '+' group operation.
	 */
	T_DERIVED
	inverse ( ) const;

	/**
	 * Inverts the element `*this` for the '+' group operation.
	 */
	void
	inverted ( )
	{ this->underlying() = this->underlying().inverse(); }

	void
	operator+= (const T_DERIVED& g);
	
	/**
	 * \return the group addition of `*this` and \p g.
	 */
	T_DERIVED
	operator+ (const T_DERIVED& g) const
	{
		T_DERIVED res(this->underlying());
		res += g;
		return res;
	}

	/**
	 * Performs `*this` + \p g.inverse().
	 */
	T_DERIVED
	operator- (const T_DERIVED& g) const
	{
		T_DERIVED res(this->underlying());
		res += g.inverse();
		return res;
	}

	/**
	 * \return the element representing the group identity for operation `+`.
	 */
	static T_DERIVED
	Zero ( );

	/* Vector space operations */
	void
	operator*= (T_SCALAR_TYPE s);
	
	/**
	 * \return the dot product of `*this` by the scalar \p s.
	 */
	T_DERIVED
	operator* (T_SCALAR_TYPE s) const
	{
		T_DERIVED res(this->underlying());
		res *= s;
		return res;
	}

	/* Lie algebra operations */

	/**
	 * Implements the non-commutative Lie bracket operation \f$[a,b]\f$ where \f$a\f$ is `*this`
	 * and \f$b\f$ is \p g.
	 * \return the bracket operation between `*this` and \p g.
	 */
	T_DERIVED
	bracket (const T_DERIVED& g) const;

	static T_DERIVED
	bracket (const T_DERIVED& g1, const T_DERIVED& g2)
	{
		return g1.bracket(g2);
	}

	/* Other operations */

	T_SCALAR_TYPE
	norm () const;

	/**
	 * Implements the exponential map \f$\exp:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 * \return the \ref Lie::SO3::Group<T> implementation of the Lie group \f$SO(3)\f$ element
	 * that represents the exponential of `*this`.
	 */
	T_GROUP
	exp ( ) const;

	NOXVector<T_DOF>
	toNOXVector ( ) const;

	static const unsigned int
	dof ()
	{
		return DOF;
	}
};

template <typename T_DERIVED, typename T_GROUP, unsigned int T_DOF, typename T_SCALAR_TYPE>
const LieAlgebraBase<T_DERIVED,T_GROUP,T_DOF,T_SCALAR_TYPE>
operator* (T_SCALAR_TYPE s, const T_DERIVED g)
{
	return g*s;
}

#endif
