#ifndef DEF_LIE_SO3_OPERATOR
#define DEF_LIE_SO3_OPERATOR

#include <Eigen/Dense>
#include "lie_SO3_group.hpp"
#include "lie_SO3_algebra.hpp"

#define CoAlgebra Algebra

namespace Lie {
namespace SO3 {

template <typename T>
class Operator {
public:
	static Group<T> AD (Group<T> const& G1, Group<T> const& G2) {
		Group<T> Res(G1*G2*G1.invert());
		return Res;
	}

	static Algebra<T> Ad (Group<T> const& G, Algebra<T> const& g) {
		Eigen::Matrix<T,3,1> V = G.matrix()*g.vect();
		Algebra<T> res(V);
		return res;
	}

	static CoAlgebra<T> Ad_star (Group<T> const& G, CoAlgebra<T> const& gs) {
		CoAlgebra<T> res(G.invert().matrix()*gs);
		return res;
	}
	
	static Algebra<T> ad (Algebra<T> const& g1, Algebra<T> const& g2) {
		return bracket(g1,g2);
	}

	static CoAlgebra<T> ad_star (Algebra<T> const& g, CoAlgebra<T> const& gs) {
		CoAlgebra<T> res(gs.vect().cross(g.vect()));
		return res;
	}
	
	static Algebra<T> bracket (Algebra<T> const& g1, Algebra<T> const& g2) {
		Algebra<T> res(g1.vect().cross(g2.vect()));
		return res;
	}

	static Group<T> cay (Algebra<T> const& g) {
		Eigen::Matrix<T,3,3> W = g.matrix();
		Eigen::Matrix<T,3,1> w = g.vect();
		Group<T> Res(Eigen::Matrix<T,3,3>::Identity()+T(4)*(W+W*W/T(2))/(T(4)+w.transpose()*w));
		Res.normalise();
		return Res;
	}

	static Algebra<T> cay_inv (Group<T> const& G) {
		// Warning : division by 0 if rotation angle is +-pi
		Eigen::Matrix<T,3,3> R(G.rotationMatrix());
		Eigen::Matrix<T,3,3> M((R-R.transpose())*(T(2)/(T(1)+R.trace())));
		Algebra<T> res(M);
		return res;
	}

	static Algebra<T> dcay (Algebra<T> const& g1, Algebra<T> const& g2) {
		Eigen::Matrix<T,3,3> W = g1.matrix();
		Eigen::Matrix<T,3,1> w = g1.vect();
		Eigen::Matrix<T,3,3> M(T(2)*(T(2)*Eigen::Matrix<T,3,3>::Identity()+W)/(T(4)+w.transpose()*w));
		Algebra<T> res(M*g2.vect());
		return res;
	}

	static Algebra<T> dcay_inv (Algebra<T> const& g1, Algebra<T> const& g2) {
		Eigen::Matrix<T,3,3> W = g1.matrix();
		Eigen::Matrix<T,3,1> w = g1.vect();
		Eigen::Matrix<T,3,3> M(Eigen::Matrix<T,3,3>::Identity()-W/T(2)+w*w.transpose()/T(4));
		Algebra<T> res(M*g2.vect());
		return res;
	}
};

} // namespace SO3
} // namespace Lie

#endif
