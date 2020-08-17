#ifndef DEF_LIE_SE3_OPERATOR
#define DEF_LIE_SE3_OPERATOR

#include <Eigen/Dense>
//#include "lie_operator_template.hpp"
#include "lie_SE3_group.hpp"
#include "lie_SE3_algebra.hpp"
#include "lie_SO3_operator.hpp"

#define CoAlgebra Algebra

namespace Lie {
namespace SE3 {

template <typename T>
class Operator {
public:
	static CoAlgebra<T> Ad_star (Group<T> const& G, CoAlgebra<T> const& gs) {
		// Ad^*(R,v)(mu,beta)=(R^T*(mu-vxbeta),R^T*beta)
		CoAlgebra<T> res(	G.rotationMatrix().transpose() \
								*(gs.rotationVect()-G.translationVect().cross(gs.translationVect())), \
							G.rotationMatrix().transpose()*gs.translationVect());
		return res;
	}

	static Group<T> cay (Algebra<T> const& g) {
		SO3::Algebra<T> om(g.getSO3());
		//Eigen::Matrix<T,3,3> M(Eigen::Matrix<T,3,3>::Identity());
		Eigen::Matrix<T,3,1> V((Eigen::Matrix<T,3,3>::Identity()+om.matrix()/T(2)+om.vect()*om.vect().transpose()/T(4)) \
									*g.translationVect()*(T(4)/(T(4)+om.vect().transpose()*om.vect())));
		Group<T> Res(SO3::Operator<T>::cay(om),V);
		return Res;
	}

	static Algebra<T> cay_inv (Group<T> G) {
		// Warning division by 0 for angle = +-pi
		Eigen::Matrix<T,3,1> V(2*(G.rotationMatrix()+Eigen::Matrix<T,3,3>::Identity()).inverse()*G.translationVect());
		Algebra<T> res(SO3::Operator<T>::cay_inv(G.getSO3()),V);
		return res;
	}

	static Algebra<T> dcay_inv_star (Algebra<T> const& g1, Algebra<T> const& g2) {
		Eigen::Matrix<T,6,6> M(Eigen::Matrix<T,6,6>::Zero());
		Eigen::Matrix<T,3,3> I(Eigen::Matrix<T,3,3>::Identity());
		SO3::Algebra<T> gam(g1.translationVect());
		M.block(0,0,3,3) = I-g1.rotationMatrix()/T(2)+g1.rotationVect()*g1.rotationVect().transpose()/T(4);
		M.block(3,0,3,3) = -(I-g1.rotationMatrix()/T(2))*gam.rotationMatrix()/T(2);
		M.block(3,3,3,3) = I-g1.rotationMatrix()/T(2);
		Algebra<T> res(M.transpose()*g2.vect());
		return res;
	}
};

} // namespace SE3
} // namespace Lie

#endif
