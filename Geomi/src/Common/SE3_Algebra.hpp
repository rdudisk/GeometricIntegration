#ifndef DEF_COMMON_SE3_ALGEBRA
#define DEF_COMMON_SE3_ALGEBRA

#include <string>

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace SE3
{

/**
 * Class Lie algebra \f$\mathfrak{se}(3)\f$ implementation.
 * 
 * \tparam T Floating point type used for internal representation of coefficients.
 */

template <typename T_SCALAR>
class Algebra : public LieAlgebraBase<	Algebra<T_SCALAR>,
										Group<T_SCALAR>,
										6,
										T_SCALAR>
{
protected:
	Eigen::Matrix<T_SCALAR,3,1> m_rot;
	Eigen::Matrix<T_SCALAR,3,1> m_trans;

public:
	Algebra<T_SCALAR> ( )
	: m_rot(Eigen::Matrix<T_SCALAR,3,1>::Zero()), m_trans(Eigen::Matrix<T_SCALAR,3,1>::Zero())
	{ }

	Algebra<T_SCALAR> (const Eigen::Matrix<T_SCALAR,3,1>& rot, const Eigen::Matrix<T_SCALAR,3,1>& trans)
	: m_rot(rot), m_trans(trans)
	{ }

	Algebra<T_SCALAR> (const Eigen::Matrix<T_SCALAR,6,1>& v)
	: m_rot(v.head(3)), m_trans(v.tail(3))
	{ }

	/* Group operations */

	/**
	 * \return the group inverse of `*this`.
	 */
	Algebra<T_SCALAR>
	inverse ( ) const
	{
		Algebra<T_SCALAR> a(-m_rot,-m_trans);
		return a;
	}

	void
	operator+= (const Algebra<T_SCALAR>& g)
	{ m_rot += g.rot(); m_trans += g.trans(); }

	/**
	 * \return the element representing the group identity for operation `+`.
	 */
	static Algebra<T_SCALAR>
	Identity ( )
	{
		Algebra<T_SCALAR> res(Eigen::Matrix<T_SCALAR,6,1>::Zero());
		return res;
	}
	
	using LieAlgebraBase<Algebra<T_SCALAR>,Group<T_SCALAR>,6,T_SCALAR>::operator+;
	using LieAlgebraBase<Algebra<T_SCALAR>,Group<T_SCALAR>,6,T_SCALAR>::operator-;

	/* Vector space operations */

	void
	operator*= (T_SCALAR s)
	{ m_rot *= s; m_trans *= s; }
	
	/**
	 * \return the dot product of `*this` by the scalar \p s.
	 */

	using LieAlgebraBase<Algebra<T_SCALAR>,Group<T_SCALAR>,6,T_SCALAR>::operator*;

	/* Lie algebra operations */

	/**
	 * Implements the non-commutative Lie bracket operation \f$[a,b]\f$ where \f$a\f$ is `*this`
	 * and \f$b\f$ is \p g.
	 * \return the bracket operation between `*this` and \p g.
	 */
	Algebra<T_SCALAR>
	bracket (const Algebra<T_SCALAR>& g) const
	{ return Algebra<T_SCALAR>(m_rot.cross(g.v()),m_rot.cross(g.trans())-g.rot().cross(m_trans)); }

	Algebra<T_SCALAR>
	Ad (const Group<T_SCALAR>& g) const
	{ return Algebra<T_SCALAR>(g.rotationMatrix()*m_rot,g.rotationMatrix()*m_trans+g.translationVector().cross(g.rotationMatrix()*m_rot)); }

	/** Implements \f$xi.Ad_star(g) = Ad^*_g(xi)\f$ */
	Algebra<T_SCALAR>
	Ad_star (const Group<T_SCALAR>& g) const
	{
		//Group<T_SCALAR> f = g.inverse();
		//return Algebra<T_SCALAR>(f.rotationMatrix()*m_rot+f.trans().cross(f.rotationMatrix()*m_trans),f.rotationMatrix()*m_trans);
		return Algebra<T_SCALAR>(g.rotationMatrix().transpose()*
				(m_rot-g.trans().cross(m_trans)),
				g.rotationMatrix().transpose()*m_trans);
	}

	/* Other operations */

	Eigen::Matrix<T_SCALAR,4,1>
	operator* (const Eigen::Matrix<T_SCALAR,4,1>& vec) const
	{ return this->toMatrix()*vec; }

	/* Accessors */

	Eigen::Matrix<T_SCALAR,3,1>
	rot ( ) const
	{ return m_rot; } 

	void
	rot (const Eigen::AngleAxis<T_SCALAR>& aa)
	{ m_rot = aa.angle()*aa.axis(); }

	void
	rot (const Eigen::Matrix<T_SCALAR,3,1>& vec)
	{ m_rot = vec; }

	Eigen::Matrix<T_SCALAR,3,1>
	trans ( ) const
	{ return m_trans; } 

	void
	trans (const Eigen::Matrix<T_SCALAR,3,1>& vec)
	{ m_trans = vec; }

	T_SCALAR const&
	operator[] (size_t index) const
	{ return (index<3) ? m_rot[index] : m_trans[index-3]; }

	T_SCALAR&
	operator[] (size_t index)
	{ return (index<3) ? m_rot[index] : m_trans[index-3]; }
	
	void
	normalizeRotation ( )
	{ m_rot.normalize(); }

	Algebra<T_SCALAR>
	normalizedRotation ( ) const
	{
		Algebra<T_SCALAR> res(*this);
		res.normalizeRotation();
		return res;
	}

	T_SCALAR
	norm ( ) const
	{
		double n1 = m_rot.norm(), n2 = m_trans.norm();
		return sqrt(n1*n1+n2*n2);
	}

	// vrai ?
	T_SCALAR
	angle ( ) const
	{ return m_rot.norm(); }

	// vrai ?
	Eigen::Matrix<T_SCALAR,3,1>
	axis ( ) const
	{ return m_rot.normalized(); }

	/**
	 * Implements the exponential map \f$\exp:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 * \return the \ref Lie::SO3::Group<T> implementation of the Lie group \f$SO(3)\f$ element
	 * that represents the exponential of `*this`.
	 */
	/*
	Group<T_SCALAR>
	exp ( ) const
	{
		T_SCALAR a = m_rot.norm()/2.0;
		Eigen::Matrix<T_SCALAR,4,1> V;
		//V << cos(a) << sin(a)*m_rot.normalized();
		//return Group<T>(V);
		return Group<T_SCALAR>(Eigen::AngleAxis<T_SCALAR>(cos(a),sin(a)*m_rot.normalized()));
	}

	Eigen::Matrix<T_SCALAR,3,3>
	partialExp (const unsigned int i) const
	{
		T_SCALAR nm = this->norm();
		Eigen::Matrix<T_SCALAR,3,3> K = (this->normalized()).toRotationMatrix();
		T_SCALAR c = cos(nm), s = sin(nm);

		Eigen::Matrix<T_SCALAR,3,3> M = Algebra<T_SCALAR>::GeneratorMatrix(i);
		Eigen::Matrix<T_SCALAR,3,3> dexpdw;

		if (isZero<T_SCALAR>(nm))
			dexpdw = M;
		else
			dexpdw = (c-(s/nm))*((*this)[i]*(1.0/nm))*K + (s/nm)*M + (s+(1.0-c)*(2.0/nm))*((*this)[i]*(1.0/nm))*K*K + ((1-c)/nm)*(M*K+K*M);

		return dexpdw;
	}

	Algebra<T_SCALAR>
	computeDExpRInv (const Algebra<T_SCALAR>& y, unsigned int order_p = 0) const
	{
		Algebra<T_SCALAR> res = y;
		if (order_p>0) {
			Algebra<T_SCALAR> A = y;
			for (int i=0; i<order_p; i++) {
				A = this->bracket(A);
				res += BERNOULLI_NUMBERS[i+1]*A;
			}
		}
		return res;
	}
	*/

	/**
	 * Implements the Cayley map \f$cay:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 * \return the \ref Lie::SO3::Group<T> implementation of the Lie group \f$SO(3)\f$ element
	 * that represents the image of `*this` by the Cayley map.
	 */
	Group<T_SCALAR>
	cay ( ) const
	{
		T_SCALAR	n =   m_rot.norm(),
						den = 4.0+n*n;
		Eigen::Matrix<T_SCALAR,3,3> W = this->rotationMatrix();
		Eigen::Matrix<T_SCALAR,4,4> M;
		M.block(0,0,3,3) = Eigen::Matrix<T_SCALAR,3,3>::Identity()
							+ (4.0/den)*(W+0.5*W*W);
		M.block(0,3,3,1) = (4.0/den)*(Eigen::Matrix<T_SCALAR,3,3>::Identity()+0.5*W+0.25*m_rot*(m_rot.transpose()))*m_trans;
		M.block(3,0,1,3) = Eigen::Matrix<T_SCALAR,1,3>::Zero();
		M(3,3) = 1.0;
		return Group<T_SCALAR>(M);
	}

	static Algebra<T_SCALAR>
	cay_inv (const Group<T_SCALAR>& g)
	{ return fromMatrix(-2.0*(Eigen::Matrix<T_SCALAR,4,4>::Identity()+g.matrix()).inverse()*(Eigen::Matrix<T_SCALAR,4,4>::Identity()-g.matrix())); }
	// Cette implementation n'a pas l'air specialement mieux
	//{ return Algebra<T_SCALAR>(SO3::Algebra<T_SCALAR>::cay_inv(SO3::Group<T_SCALAR>(g.q())).toVector(), 2.0*(g.rotationMatrix()+Eigen::Matrix<T_SCALAR,3,3>::Identity()).inverse()*g.trans()); }

	Eigen::Matrix<T_SCALAR,6,6>
	dCayRInv () const
	{
		Eigen::Matrix<T_SCALAR,6,6> M = Eigen::Matrix<T_SCALAR,6,6>::Identity();
		M.block(0,0,3,3) += -0.5*this->rotationMatrix()+0.25*m_rot*(m_rot).transpose();
		M.block(3,0,3,3) = -0.5*(Eigen::Matrix<T_SCALAR,3,3>::Identity()-0.5*this->rotationMatrix())*toRotationMatrix(m_trans);
		M.block(3,3,3,3) -= 0.5*this->rotationMatrix();
		return M;
		/* Autre implémentation
		Eigen::Matrix<T_SCALAR,6,6> M = Eigen::Matrix<T_SCALAR,6,6>::Zero();
		Eigen::Matrix<T_SCALAR,3,3> I = Eigen::Matrix<T_SCALAR,3,3>::Identity();
		M.block(0,0,3,3) = I-0.5*this->rotationMatrix()+0.25*m_rot*(m_rot.transpose());
		M.block(3,0,3,3) = -0.5*(I-0.5*this->rotationMatrix())*toRotationMatrix(m_trans);
		M.block(3,3,3,3) = I-0.5*this->rotationMatrix();
		return M;*/
	}

	Algebra<T_SCALAR>
	dCayRInv (const Algebra<T_SCALAR>& g) const
	{ return Algebra<T_SCALAR>(this->dCayRInv()*g.toVector()); }

	Eigen::Matrix<T_SCALAR,6,6>
	dCayRInvStar () const
	{ return (this->dCayRInv()).transpose(); }

	Algebra<T_SCALAR>
	dCayRInvStar (const Algebra<T_SCALAR>& g) const
	{ return Algebra<T_SCALAR>(this->dCayRInvStar()*g.toVector()); }

	/**
	 * \return the axis-angle representation of `*this`.
	 */
	/*
	Eigen::AngleAxis<T_SCALAR>
	toAngleAxis ( ) const
	{ return Eigen::AngleAxis<T_SCALAR>(m_rot.norm(),m_rot.normalized()); }
	*/

	Eigen::Matrix<T_SCALAR,3,3>
	toRotationMatrix (const Eigen::Matrix<T_SCALAR,3,1>& v) const
	{
		Eigen::Matrix<T_SCALAR,3,3> mat;
		mat <<	T_SCALAR(0),		-v[2],	v[1],
				v[2],		T_SCALAR(0),		-v[0],
				-v[1],	v[0],		T_SCALAR(0);
		return mat;
	}

	Eigen::Matrix<T_SCALAR,3,3>
	rotationMatrix ( ) const
	{ return this->toRotationMatrix(this->m_rot); }

	static Algebra<T_SCALAR>
	fromMatrix (const Eigen::Matrix<T_SCALAR,4,4>& m)
	{ return Algebra<T_SCALAR>(Eigen::Matrix<T_SCALAR,3,1>(m(2,1),m(0,2),m(1,0)),m.block(0,3,3,1)); }

	Eigen::Matrix<T_SCALAR,6,1>
	toVector ( ) const
	{ 
		Eigen::Matrix<T_SCALAR,6,1> ret;
		ret.head(3) = m_rot;
		ret.tail(3) = m_trans;
		return ret;
	}

	NOXVector<6>
	toNOXVector ( ) const
	{ return NOXVector<6>(this->toVector()); }

	Eigen::Matrix<T_SCALAR,4,4>
	toMatrix ( ) const
	{
		Eigen::Matrix<T_SCALAR,4,4> M = Eigen::Matrix<T_SCALAR,4,4>::Zero();
		M.block(0,0,3,3) = this->rotationMatrix();
		M.block(0,3,3,1) = this->m_trans;
		return M;
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
	static Eigen::Matrix<T_SCALAR,4,4>
	GeneratorMatrix (int const i)
	{
		Eigen::Matrix<T_SCALAR,4,4> res(Eigen::Matrix<T_SCALAR,4,4>::Zero());
		switch (i) {
			case 0:
				res(1,2) = T_SCALAR(-1); res(2,1) = T_SCALAR(1);
				break;
			case 1:
				res(0,2) = T_SCALAR(1);  res(2,0) = T_SCALAR(-1);
				break;
			case 2:
				res(0,1) = T_SCALAR(-1); res(1,0) = T_SCALAR(1);
				break;
			case 3:
				res(0,3) = T_SCALAR(1);
				break;
			case 4:
				res(1,3) = T_SCALAR(1);
				break;
			case 5:
				res(2,3) = T_SCALAR(1);
				break;
			default:
				break;
		}
		return res;
	}

	static Eigen::Matrix<T_SCALAR,6,1>
	GeneratorVector (int const i)
	{
		Eigen::Matrix<T_SCALAR,6,1> res(Eigen::Matrix<T_SCALAR,6,1>::Zero());
		res(i) = T_SCALAR(1);
		return res;
	}

	static Algebra<T_SCALAR>
	Generator (int const i)
	{ return Algebra(Algebra::GeneratorVector(i)); }

	static Algebra<T_SCALAR>
	Zero ( )
	{ return Algebra(Eigen::Matrix<T_SCALAR,6,1>::Zero()); }

	static Algebra<T_SCALAR>
	Random ( )
	{ return Algebra(Eigen::Matrix<T_SCALAR,6,1>::Random()); }
};

} // namespace SO3

template <typename T_SCALAR>
std::ostream&
operator<< (std::ostream& stream, SE3::Algebra<T_SCALAR> const& g)
{
	stream << g.toVector();
	return stream;
}

template <typename T_SCALAR>
const SE3::Algebra<T_SCALAR>
operator* (T_SCALAR s, const SE3::Algebra<T_SCALAR> g)
{ return g*s; }

#endif
