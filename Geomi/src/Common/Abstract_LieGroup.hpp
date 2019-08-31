#ifndef DEF_COMMON_ABSTRACT_LIEGROUP
#define DEF_COMMON_ABSTRACT_LIEGROUP

namespace Abstract
{

class LieGroup {
public:
	static const unsigned int DOF;

public:
	LieGroup ( )
	{ }

	virtual 
	~LieGroup ( )
	{ }

	/**
	 * @return the group inverse of `*this`.
	 */
	virtual LieGroup&
	inverse ( ) const = 0;

	/**
	 * Inverts the group element `*this`
	 */
	void
	inverted ( )
	{ *this=this->inverse(); }

	/**
	 * \return the element representing the group identity for operation `*`.
	 */
	static LieGroup&
	Identity ( );

	/**
	 * Group operation '*'.
	 */
	virtual void
	operator*= (LieGroup const& g) = 0;

	/**
	 * \return the vector representation of the rotation, that is the vector \f$\vec v\f$
	 * such that the rotation of any given vector \f$\vec u\f$ is the result of \f$vec v\wedge\vec u\f$.
	 */
	virtual NOXVector<DOF>
	toNOXVector ( ) const = 0;

	static const unsigned int
	dof ()
	{
		return DOF;
	}
};

} // namespace Abstract

#endif
