#ifndef DEF_COMMON_SO3_ALGEBRA
#define DEF_COMMON_SO3_ALGEBRA

#include <string>

#include <Eigen/Dense>
#include <Eigen/Geometry>

//#include "Common_SO3_Group.hpp"
//#include "include/Common/NOXVector.hpp"
//#include "include/Common/Utils.hpp"

namespace SO3
{

/**
 * Class Lie algebra \f$\mathfrak{so}(3)\f$ implementation.
 * \tparam T Floating point type used for internal representation of coefficients.
 */

template <typename T_SCALAR_TYPE>
class Algebra : public LieAlgebraBase<	Algebra<T_SCALAR_TYPE>,
										Group<T_SCALAR_TYPE>,
										3,
										T_SCALAR_TYPE>
{
protected:
	Eigen::Matrix<T_SCALAR_TYPE,3,1> m_v;

public:
	Algebra<T_SCALAR_TYPE> ( )
	: m_v(Eigen::Matrix<T_SCALAR_TYPE,3,1>::Zero())
	{ }

	Algebra<T_SCALAR_TYPE> (const Eigen::Matrix<T_SCALAR_TYPE,3,1>& v)
	: m_v(v)
	{ }

	Algebra<T_SCALAR_TYPE> (const T_SCALAR_TYPE a, const T_SCALAR_TYPE b, const T_SCALAR_TYPE c)
	: m_v(Eigen::Matrix<T_SCALAR_TYPE,3,1>(a,b,c))
	{ }

	// Removed because ambiguity with  Algebra<T_SCALAR_TYPE> (const Eigen::Matrix<T_SCALAR_TYPE,3,1>& v)
	/* Algebra<T_SCALAR_TYPE> (const NOXVector<3> v)
	: m_v(v)
	{ }
	*/

	/* Group operations */

	/**
	 * \return the group inverse of `*this`.
	 */
	Algebra<T_SCALAR_TYPE>
	inverse ( ) const
	{
		Algebra<T_SCALAR_TYPE> a(-m_v);
		return a;
	}

	void
	operator+= (const Algebra<T_SCALAR_TYPE>& g)
	{ m_v += g.v(); }

	/**
	 * \return the element representing the group identity for operation `+`.
	 */
	static Algebra<T_SCALAR_TYPE>
	Identity ( )
	{
		Algebra<T_SCALAR_TYPE> res(Eigen::Matrix<T_SCALAR_TYPE,3,1>::Zero());
		return res;
	}
	
	using LieAlgebraBase<Algebra<T_SCALAR_TYPE>,Group<T_SCALAR_TYPE>,3,T_SCALAR_TYPE>::operator+;
	using LieAlgebraBase<Algebra<T_SCALAR_TYPE>,Group<T_SCALAR_TYPE>,3,T_SCALAR_TYPE>::operator-;

	/* Vector space operations */

	void
	operator*= (T_SCALAR_TYPE s)
	{ m_v *= s; }
	
	/**
	 * \return the dot product of `*this` by the scalar \p s.
	 */

	using LieAlgebraBase<Algebra<T_SCALAR_TYPE>,Group<T_SCALAR_TYPE>,3,T_SCALAR_TYPE>::operator*;
	/*
	Algebra<T_SCALAR_TYPE>
	operator* (const T_SCALAR_TYPE& s) const
	{
		return LieAlgebraBase<Algebra<T_SCALAR_TYPE>,Group<T_SCALAR_TYPE>,3,T_SCALAR_TYPE>::operator*(s);
	}
	*/

	/* Lie algebra operations */

	/**
	 * Implements the non-commutative Lie bracket operation \f$[a,b]\f$ where \f$a\f$ is `*this`
	 * and \f$b\f$ is \p g.
	 * \return the bracket operation between `*this` and \p g.
	 */
	Algebra<T_SCALAR_TYPE>
	bracket (const Algebra<T_SCALAR_TYPE>& g) const
	{ return this->m_v.cross(g.v()); }

	Algebra<T_SCALAR_TYPE>
	Ad (const Group<T_SCALAR_TYPE>& g) const
	{ return Algebra<T_SCALAR_TYPE>(g.toRotationMatrix()*this->toVector()); }

	Algebra<T_SCALAR_TYPE>
	Ad_star (const Group<T_SCALAR_TYPE>& g) const
	{ return Algebra<T_SCALAR_TYPE>(g.toRotationMatrix().transpose()*this->toVector()); }

	/* Other operations */

	/**
	 * Transforms the size 3 vector \p vec by the algebra element represented by `*this`.
	 * With a skew matrix representation, it implements the operation \f$\widehat\omega v\f$;
	 * with a vector reprentation it implements \f$\omega\wedge v\f$.
	 * \return the result of the action of `*this` on \p vec.
	 */
	Eigen::Matrix<T_SCALAR_TYPE,3,1>
	operator* (const Eigen::Matrix<T_SCALAR_TYPE,3,1>& vec) const
	{ return m_v.cross(vec); }

	/* Accessors */

	Eigen::Matrix<T_SCALAR_TYPE,3,1>
	v ( ) const
	{ return m_v; } 

	void
	v (const Eigen::AngleAxis<T_SCALAR_TYPE>& aa)
	{ m_v = aa.angle()*aa.axis(); }

	void
	v (const Eigen::Matrix<T_SCALAR_TYPE,3,1>& vec)
	{ m_v = vec; }

	T_SCALAR_TYPE const&
	operator[] (size_t index) const
	{ return m_v[index]; }

	T_SCALAR_TYPE&
	operator[] (size_t index)
	{ return m_v[index]; }
	
	// Les fonctions de normalisation sont-elles vraiment utiles ?
	void
	normalize ( )
	{ m_v.normalize(); }

	Algebra<T_SCALAR_TYPE>
	normalized ( ) const
	{
		Algebra<T_SCALAR_TYPE> res(*this);
		res.normalize();
		return res;
	}

	T_SCALAR_TYPE
	norm ( ) const
	{ return m_v.norm(); }

	T_SCALAR_TYPE
	angle ( ) const
	{ return m_v.norm(); }

	Eigen::Matrix<T_SCALAR_TYPE,3,1>
	axis ( ) const
	{ return m_v.normalized(); }

	/**
	 * Implements the exponential map \f$\exp:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 * \return the \ref Lie::SO3::Group<T> implementation of the Lie group \f$SO(3)\f$ element
	 * that represents the exponential of `*this`.
	 */
	Group<T_SCALAR_TYPE>
	exp ( ) const
	{
		T_SCALAR_TYPE a = m_v.norm()/2.0;
		Eigen::Matrix<T_SCALAR_TYPE,4,1> V;
		//V << cos(a) << sin(a)*m_v.normalized();
		//return Group<T>(V);
		return Group<T_SCALAR_TYPE>(Eigen::AngleAxis<T_SCALAR_TYPE>(cos(a),sin(a)*m_v.normalized()));
	}

	Eigen::Matrix<T_SCALAR_TYPE,3,3>
	partialExp (const unsigned int i) const
	{
		T_SCALAR_TYPE nm = this->norm();
		Eigen::Matrix<T_SCALAR_TYPE,3,3> K = (this->normalized()).toRotationMatrix();
		T_SCALAR_TYPE c = cos(nm), s = sin(nm);

		Eigen::Matrix<T_SCALAR_TYPE,3,3> M = Algebra<T_SCALAR_TYPE>::GeneratorMatrix(i);
		Eigen::Matrix<T_SCALAR_TYPE,3,3> dexpdw;

		if (isZero<T_SCALAR_TYPE>(nm))
			dexpdw = M;
		else
			dexpdw = (c-(s/nm))*((*this)[i]*(1.0/nm))*K + (s/nm)*M + (s+(1.0-c)*(2.0/nm))*((*this)[i]*(1.0/nm))*K*K + ((1-c)/nm)*(M*K+K*M);

		return dexpdw;
	}

	Algebra<T_SCALAR_TYPE>
	computeDExpRInv (const Algebra<T_SCALAR_TYPE>& y, unsigned int order_p = 0) const
	{
		Algebra<T_SCALAR_TYPE> res = y;
		if (order_p>0) {
			Algebra<T_SCALAR_TYPE> A = y;
			for (int i=0; i<order_p; i++) {
				A = this->bracket(A);
				res += BERNOULLI_NUMBERS[i+1]*A;
			}
		}
		return res;
	}

	/**
	 * Implements the Cayley map \f$cay:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 * \return the \ref Lie::SO3::Group<T> implementation of the Lie group \f$SO(3)\f$ element
	 * that represents the image of `*this` by the Cayley map.
	 */
	Group<T_SCALAR_TYPE>
	cay ( ) const
	{
		T_SCALAR_TYPE	n =   m_v.norm(),
						den = 4.0+n*n;
		Eigen::Matrix<T_SCALAR_TYPE,4,1> V;
		//V << 1.0 - 2.0*n*n/den << m_v.normalized() * 4.0*n/den;
		//return Group<T>(V);
		return Group<T_SCALAR_TYPE>(Eigen::AngleAxis<T_SCALAR_TYPE>(1.0-2.0*n*n/den,m_v.normalized()*4.0*n/den));
	}

	static Algebra<T_SCALAR_TYPE>
	cay_inv (const Group<T_SCALAR_TYPE>& g)
	{ return fromRotationMatrix(2.0*(g.toRotationMatrix()-Eigen::Matrix<T_SCALAR_TYPE,3,3>::Identity())*(g.toRotationMatrix()+Eigen::Matrix<T_SCALAR_TYPE,3,3>::Identity()).inverse()); }

	Eigen::Matrix<T_SCALAR_TYPE,3,3>
	dCayRInv () const
	{ return Eigen::Matrix<T_SCALAR_TYPE,3,3>::Identity()-0.5*this->toRotationMatrix()+0.25*this->toVector()*(this->toVector().transpose()); }

	Algebra<T_SCALAR_TYPE>
	dCayRInv (const Algebra<T_SCALAR_TYPE>& g) const
	{ return Algebra<T_SCALAR_TYPE>(this->dCayRInv()*g.toVector()); }

	/**
	 * \return the axis-angle representation of `*this`.
	 */
	Eigen::AngleAxis<T_SCALAR_TYPE>
	toAngleAxis ( ) const
	{ return Eigen::AngleAxis<T_SCALAR_TYPE>(m_v.norm(),m_v.normalized()); }

	Eigen::Matrix<T_SCALAR_TYPE,3,3>
	toRotationMatrix ( ) const
	{
		Eigen::Matrix<T_SCALAR_TYPE,3,3> mat;
		mat <<	T_SCALAR_TYPE(0),		-m_v[2],	m_v[1],
				m_v[2],		T_SCALAR_TYPE(0),		-m_v[0],
				-m_v[1],	m_v[0],		T_SCALAR_TYPE(0);
		return mat;
	}

	static Algebra<T_SCALAR_TYPE>
	fromRotationMatrix (const Eigen::Matrix<T_SCALAR_TYPE,3,3>& m)
	{
		return Algebra<T_SCALAR_TYPE>(m(2,1),m(0,2),m(1,0));
	}

	Eigen::Matrix<T_SCALAR_TYPE,3,1>
	toVector ( ) const
	{ return m_v; }

	NOXVector<3>
	toNOXVector ( ) const
	{ return NOXVector<3>(m_v); }

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
	static Eigen::Matrix<T_SCALAR_TYPE,3,3>
	GeneratorMatrix (int const i)
	{
		Eigen::Matrix<T_SCALAR_TYPE,3,3> res(Eigen::Matrix<T_SCALAR_TYPE,3,3>::Zero());
		if (i == 0)			{ res(1,2) = T_SCALAR_TYPE(-1); res(2,1) = T_SCALAR_TYPE(1);  }
		else if (i == 1)	{ res(0,2) = T_SCALAR_TYPE(1);  res(2,0) = T_SCALAR_TYPE(-1); }
		else				{ res(0,1) = T_SCALAR_TYPE(-1); res(1,0) = T_SCALAR_TYPE(1);  }
		return res;
	}

	static Eigen::Matrix<T_SCALAR_TYPE,3,1>
	GeneratorVector (int const i)
	{
		Eigen::Matrix<T_SCALAR_TYPE,3,1> res(Eigen::Matrix<T_SCALAR_TYPE,3,1>::Zero());
		res(i) = T_SCALAR_TYPE(1);
		return res;
	}

	static Algebra<T_SCALAR_TYPE>
	Generator (int const i)
	{ return Algebra(Algebra::GeneratorVector(i)); }

	static Algebra<T_SCALAR_TYPE>
	Zero ( )
	{ return Algebra(Eigen::Matrix<T_SCALAR_TYPE,3,1>::Zero()); }
};

} // namespace SO3

template <typename T_SCALAR_TYPE>
std::ostream&
operator<< (std::ostream& stream, SO3::Algebra<T_SCALAR_TYPE> const& g)
{
	stream << g.v();
	return stream;
}

template <typename T_SCALAR_TYPE>
const SO3::Algebra<T_SCALAR_TYPE>
operator* (T_SCALAR_TYPE s, const SO3::Algebra<T_SCALAR_TYPE> g)
{ return g*s; }

#endif
