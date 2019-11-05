#ifndef DEF_COMMON_LIEALGEBRABASE
#define DEF_COMMON_LIEALGEBRABASE

/**
 * This class defines a Lie algebra \f$(\mathfrak g,[\cdot,\cdot])\f$
 * where \f$(\mathfrak g,\mathbb F,+,\cdot)\f$ is a vector space
 * over the field \f$\mathbb F\f$,
 * with the vector space addition
 * \f$ +:\mathfrak g\times\mathfrak g\rightarrow\mathfrak g\f$.
 * and the scalar multiplication
 * \f$\cdot:\mathbb F\times\mathfrak g\rightarrow\mathfrak g\f$.
 * Here the field \f$\mathbb F\f$ is represented by the templated type
 * `T_SCALAR_TYPE`.
 *
 * The design of the class is based on the Curious Recurring Template Pattern.
 * The class is intended to be inherited.
 * The function members to be implemented by the user are the following:
 *		- `operator+= (const T_DERIVED&)`
 *		- `inverse ()`
 *		- `Zero ()`
 *		- `operator*= (T_SCALAR_TYPE)`
 *		- `bracket (const T_DERIVED&g)`
 *		- `Ad (const T_GROUP&)`
 *		- `Ad_star (const T_GROUP&)`
 *		- `norm ()`
 *		- `toNOXVector ()`
 */
template <typename T_DERIVED, typename T_GROUP, unsigned int T_DOF, typename T_SCALAR_TYPE>
class LieAlgebraBase : public CRTP<T_DERIVED>
{
public:
	static const unsigned int DOF = T_DOF;

public:
	/* Group operations */

	/**
	 * This function has to be implemented by the user.
	 * It should return the mathematical inverse of `*this`
	 * with respect to the vector space addition.
	 */
	T_DERIVED
	inverse ( ) const;

	/**
	 * Inverts `*this` in place.
	 * \see inverse()
	 */
	void
	inverted ( )
	{ this->underlying() = this->underlying().inverse(); }

	/**
	 * This function has to be implemented by the user.
	 * It should perform the vector space addition of
	 * the argument \p g and `*this`.
	 * The operation must be done in place.
	 */
	void
	operator+= (const T_DERIVED& g);
	
	/**
	 * Performs the vector space addition of `*this` and \p g.
	 * \see operator+=(const T_DRIVED& g)
	 */
	T_DERIVED
	operator+ (const T_DERIVED& g) const
	{
		T_DERIVED res(this->underlying());
		res += g;
		return res;
	}

	/**
	 * Adds `*this` to the vector space addition inverse of \p g.
	 */
	T_DERIVED
	operator- (const T_DERIVED& g) const
	{
		T_DERIVED res(this->underlying());
		res += g.inverse();
		return res;
	}

	/**
	 * This function has to be implemented by the user.
	 * It should return the Lie algebra representing the identity element for
	 * the vector space addition.
	 */
	static T_DERIVED
	Zero ( );

	/**
	 * This function has to be implemented by the user.
	 * It should perform the scalar multiplication
	 * of the argument \p s and `*this`.
	 * The operation must be done in place.
	 */
	void
	operator*= (T_SCALAR_TYPE s);
	
	/**
	 * Performs the scalar multiplication of \p g and `*this`.
	 * \see operator*=(T_SCALAR_TYPE s)
	 */
	T_DERIVED
	operator* (const T_SCALAR_TYPE& s) const
	{
		T_DERIVED res(this->underlying());
		res *= s;
		return res;
	}

	/**
	 * This function has to be implemented by the user.
	 * It should perform the Lie bracket operation
	 * \f$[a,b]\f$ where \f$a\f$ is `*this` and \f$b\f$ is \p g.
	 */
	T_DERIVED
	bracket (const T_DERIVED& g) const;

	/**
	 * Static version of `bracket(const T_DERIVED&)`.
	 * \see bracket(const T_DERIVED&)
	 */
	// TODO: static_ added because the compiler doesn't seem to see this method
	// when T_DERIVED already exists
	static T_DERIVED
	static_bracket (const T_DERIVED& g1, const T_DERIVED& g2)
	{ return g1.bracket(g2); }

	/* I didn't go further for the documentation */

	T_DERIVED
	Ad (const T_GROUP& g) const;

	static T_DERIVED
	static_Ad (const T_DERIVED& a, const T_GROUP& g)
	{ return a.Ad(g); }

	T_DERIVED
	Ad_star (const T_GROUP& g) const;

	static T_DERIVED
	static_Ad_star (const T_GROUP& g, const T_DERIVED& a)
	{ return a.Ad_star(g); }

	/* Other operations */

	T_SCALAR_TYPE
	norm () const;

	/**
	 * Implements the exponential map \f$\exp:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 */
	T_GROUP
	exp ( ) const;

	static T_GROUP
	exp (const T_DERIVED& g)
	{ return g.exp(); }

	NOXVector<T_DOF>
	toNOXVector ( ) const;

	static const unsigned int
	dof ()
	{ return DOF; }
};

/*
 * Faire quelque chose à propos de ça (friend?)
 *
template <typename T_DERIVED, typename T_SCALAR_TYPE>
const T_DERIVED
operator* (T_SCALAR_TYPE s, const T_DERIVED g)
{
	return g*s;
}
*/

#endif
