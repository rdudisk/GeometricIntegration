#ifndef DEF_COMMON_SE3_CAY
#define DEF_COMMON_SE3_CAY

#include <Eigen/Dense>

namespace SE3
{
namespace Cay
{
	/**
	 * Implements the Cayley map \f$cay:\mathfrak{so}(3)\rightarrow SO(3)\f$ evalutated at `*this`.
	 * \return the \ref Lie::SO3::Group<T> implementation of the Lie group \f$SO(3)\f$ element
	 * that represents the image of `*this` by the Cayley map.
	 */
	// tested (used in program and seems to work)
	template <typename T_SCALAR>
	SE3::Group<T_SCALAR>
	eval (const SE3::Algebra<T_SCALAR>& a)
	{
		T_SCALAR n = a.rotationVector().norm(), den = 4.0+n*n;
		Eigen::Matrix<T_SCALAR,3,3> W = a.rotationMatrix();
		Eigen::Matrix<T_SCALAR,4,4> M = Eigen::Matrix<T_SCALAR,4,4>::Identity();
		M.block(0,0,3,3) = Eigen::Matrix<T_SCALAR,3,3>::Identity() + (4.0/den)*(W+0.5*W*W);
		M.block(0,3,3,1) = (4.0/den)*(Eigen::Matrix<T_SCALAR,3,3>::Identity()+0.5*W+0.25*a.rotationVector()*(a.rotationVector().transpose()))*a.translationVector();
		return Group<T_SCALAR>(M);
	}

	// tested (used in program and seems to work)
	template <typename T_SCALAR>
	SE3::Algebra<T_SCALAR>
	inv (const SE3::Group<T_SCALAR>& g)
	{ return SE3::Algebra<T_SCALAR>::fromMatrix(-2.0*(Eigen::Matrix<T_SCALAR,4,4>::Identity()+g.matrix()).inverse()*(Eigen::Matrix<T_SCALAR,4,4>::Identity()-g.matrix())); }
	// Cette implementation n'a pas l'air specialement mieux
	//{ return Algebra<T_SCALAR>(SO3::Algebra<T_SCALAR>::cay_inv(SO3::Group<T_SCALAR>(g.q())).toVector(), 2.0*(g.rotationMatrix()+Eigen::Matrix<T_SCALAR,3,3>::Identity()).inverse()*g.trans()); }

	template <typename T_SCALAR>
	Eigen::Matrix<T_SCALAR,6,6>
	inv_right_diff (const SE3::Algebra<T_SCALAR>& a)
	{
		Eigen::Matrix<T_SCALAR,6,6> M = Eigen::Matrix<T_SCALAR,6,6>::Identity();
		Eigen::Matrix<T_SCALAR,3,3> W = a.rotationMatrix();
		M.block(0,0,3,3) += -0.5*W+0.25*a.rotationVector()*(a.rotationVector().transpose());
		M.block(3,0,3,3) = -0.5*(Eigen::Matrix<T_SCALAR,3,3>::Identity()-0.5*W)*SE3::Algebra<T_SCALAR>::toRotationMatrix(a.translationVector());
		M.block(3,3,3,3) -= 0.5*W;
		return M;
		/* Autre implémentation
		Eigen::Matrix<T_SCALAR,6,6> M = Eigen::Matrix<T_SCALAR,6,6>::Zero();
		Eigen::Matrix<T_SCALAR,3,3> I = Eigen::Matrix<T_SCALAR,3,3>::Identity();
		M.block(0,0,3,3) = I-0.5*this->rotationMatrix()+0.25*m_rot*(m_rot.transpose());
		M.block(3,0,3,3) = -0.5*(I-0.5*this->rotationMatrix())*toRotationMatrix(m_trans);
		M.block(3,3,3,3) = I-0.5*this->rotationMatrix();
		return M;*/
	}

	template <typename T_SCALAR>
	SE3::Algebra<T_SCALAR>
	inv_right_diff (const SE3::Algebra<T_SCALAR>& a, const SE3::Algebra<T_SCALAR>& b)
	{ return SE3::Algebra<T_SCALAR>(inv_right_diff<T_SCALAR>(a)*b.vector()); }

	template <typename T_SCALAR>
	Eigen::Matrix<T_SCALAR,6,6>
	inv_right_diff_star (const SE3::Algebra<T_SCALAR>& a)
	{ return (inv_right_diff<T_SCALAR>(a)).transpose(); }

	// tested (used in program and seems to work)
	template <typename T_SCALAR>
	SE3::Algebra<T_SCALAR>
	inv_right_diff_star (const SE3::Algebra<T_SCALAR>& a, const SE3::Algebra<T_SCALAR>& b)
	{ return SE3::Algebra<T_SCALAR>(inv_right_diff_star<T_SCALAR>(a)*b.vector()); }

} // namespace Cay
} // namespace SE3

#endif
