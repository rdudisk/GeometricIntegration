#ifndef DEF_COMMON_LIEGROUPBASE
#define DEF_COMMON_LIEGROUPBASE

/**
 * This class uses the Curiously Recurring Template Pattern (CRTP).
 * To define an inherited class, do something like
 * class Derived : LieGroupBase<Derived,DOF>
 */
template <typename T_DERIVED, unsigned int T_DOF>
class LieGroupBase : CRTP<T_DERIVED>
{
public:
	static const unsigned int DOF = T_DOF;

public:
	/*LieGroupBase<T_DERIVED,T_DOF> ( )
	{ }

	virtual 
	~LieGroupBase<T_DERIVED,T_DOF> ( )
	{ }*/

	/**
	 * @return the group inverse of `*this`.
	 */
	T_DERIVED
	inverse ( ) const;

	/**
	 * Inverts the group element `*this`
	 */
	void
	inverted ( )
	{ this->underlying() = this->underlying().inverse(); }

	/**
	 * \return the element representing the group identity for operation `*`.
	 */
	static T_DERIVED
	Identity ( );

	/**
	 * Group operation '*'.
	 */
	void
	operator*= (T_DERIVED const& g);

	/**
	 * \return the vector representation of the rotation, that is the vector \f$\vec v\f$
	 * such that the rotation of any given vector \f$\vec u\f$ is the result of \f$vec v\wedge\vec u\f$.
	 */
	NOXVector<T_DOF>
	toNOXVector ( ) const;

	static const unsigned int
	dof ()
	{
		return DOF;
	}
};

#endif
