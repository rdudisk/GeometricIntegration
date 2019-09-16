#ifndef DEF_COMMON_NOXVECTOR
#define DEF_COMMON_NOXVECTOR

#include <cmath>

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Random.H"

#include <Eigen/Dense>
#include <Eigen/Geometry>

template <unsigned int T_DOF>
class NOXVector :	public virtual NOX::Abstract::Vector,
					public Eigen::Matrix<double,T_DOF,1>
{
public:
	static const unsigned int DOF = T_DOF;

public:
	~NOXVector<T_DOF> ()
	{ }

	/* See Eigen doc on why those functions exist */
	NOXVector<T_DOF> ()
	: Eigen::Matrix<double,T_DOF,1>()
	{ }

	template <typename OtherDerived>
	NOXVector<T_DOF> (const Eigen::MatrixBase<OtherDerived>& other)
	: Eigen::Matrix<double,T_DOF,1>(other)
	{ }

	template <typename OtherDerived>
	NOX::Abstract::Vector&
	operator= (const Eigen::MatrixBase<OtherDerived>& other)
	{ this->Eigen::Matrix<double,T_DOF,1>::operator=(other); return *this; }
	/* End 'See Eigen' */

	template <typename OtherDerived>
	NOX::Abstract::Vector&
	operator+= (const Eigen::MatrixBase<OtherDerived>& other)
	{ this->Eigen::Matrix<double,T_DOF,1>::operator+=(other); return *this; }

	// Copy constructor with NOX::CopyType : just ignore type
	NOXVector<T_DOF> (const NOXVector& source, NOX::CopyType type)
	: NOXVector<T_DOF>(source)
	{ }

	NOX::Abstract::Vector&
	abs (const NOXVector<T_DOF>& y)
	{ *this = y.cwiseAbs(); return *this; }

	NOX::Abstract::Vector&
	abs (const NOX::Abstract::Vector& y)
	{ return abs(dynamic_cast<const NOXVector<T_DOF>&>(y)); }

	Teuchos::RCP<NOX::Abstract::Vector>
	clone (NOX::CopyType type=NOX::DeepCopy) const
	{
		Teuchos::RCP<NOX::Abstract::Vector> tmp;
		tmp = Teuchos::rcp(new NOXVector<T_DOF>(*this,type));
		return tmp;
	}

	// TODO: (mais pas necessaire pour faire tourner la machine)
	// override createMultiVector()
	// pour ca il faut une implementation de NOX::Abstract::MultiVector, qui par defaut (le cas ici) est NOX::MultiVector

	double
	innerProduct (const NOXVector<T_DOF>& y) const
	{ return this->dot(y); }

	double
	innerProduct (const NOX::Abstract::Vector& y) const
	{ return innerProduct(dynamic_cast<const NOXVector<T_DOF>&>(y)); }

	NOX::Abstract::Vector&
	init (double gamma)
	{ this->setConstant(gamma); return *this; }

	NOX::size_type
	length () const
	{ return T_DOF; }

	double
	norm (NOX::Abstract::Vector::NormType type=NOX::Abstract::Vector::TwoNorm) const
	{
		double result;

		switch (type) {
			case NOX::Abstract::Vector::MaxNorm:
				result = this->Eigen::Matrix<double,T_DOF,1>::maxCoeff();
				break;
			case NOX::Abstract::Vector::OneNorm:
				// TODO norme 1
				//result = this.lpNorm<1>();
				//break;
			case NOX::Abstract::Vector::TwoNorm:
				result = this->Eigen::Matrix<double,T_DOF,1>::norm();
				break;
			default:
				result = this->Eigen::Matrix<double,T_DOF,1>::norm();
				break;
		}
		
		return result;
	}

	double
	norm (const NOXVector& weights) const
	{ return sqrt(((this->cwiseProduct(*this)).cwiseProduct(weights)).sum()); }

	double
	norm (const NOX::Abstract::Vector& weights) const
	{ return norm(dynamic_cast<const NOXVector<T_DOF>&>(weights)); }

	void
	print (std::ostream& stream) const
	{
		//TODO
		return;
	}

	NOX::Abstract::Vector&
	random (bool useSeed=false, int seed=1)
	{
		if (useSeed)
			NOX::Random::setSeed(seed);
		
		for (int i=0; i<T_DOF; i++)
			(*this)[i] = NOX::Random::number();

		return *this;
	}

	// TODO: definir operator=(const NOXVector&), ou bien est-ce que c'est directement herite de Eigen ?
	NOXVector<T_DOF>&
	operator= (const NOXVector<T_DOF>& y)
	{ this->Eigen::Matrix<double,T_DOF,1>::operator=(y); return *this; }

	NOX::Abstract::Vector&
	operator= (const NOX::Abstract::Vector& y)
	{ return operator=(dynamic_cast<const NOXVector<T_DOF>&>(y)); }

	NOX::Abstract::Vector&
	reciprocal (const NOXVector& y)
	{ *this = this->cwiseInverse(); return *this; }

	NOX::Abstract::Vector&
	reciprocal (const NOX::Abstract::Vector& y)
	{ return reciprocal(dynamic_cast<const NOXVector<T_DOF>&>(y)); }

	NOX::Abstract::Vector&
	scale (double gamma)
	{ *this = (*this)*gamma; return *this; }

	NOX::Abstract::Vector&
	scale (const NOXVector<T_DOF>& a)
	{ *this = this->cwiseProduct(a); return *this; }

	NOX::Abstract::Vector&
	scale (const NOX::Abstract::Vector& a)
	{ return scale(dynamic_cast<const NOXVector<T_DOF>&>(a)); }

	NOX::Abstract::Vector&
	update (double alpha, const NOXVector<T_DOF>& a, double gamma=0.0)
	{ *this = alpha*a + (*this)*gamma; return *this; }

	NOX::Abstract::Vector&
	update (double alpha, const NOX::Abstract::Vector& a, double gamma=0.0)
	{ return update(alpha,dynamic_cast<const NOXVector<T_DOF>&>(a),gamma); }

	NOX::Abstract::Vector&
	update (	double alpha,
				const NOXVector<T_DOF>& a,
				double beta,
				const NOXVector<T_DOF>& b,
				double gamma=0.0)
	{ *this = alpha*a + beta*b + (*this)*gamma; return *this; }

	NOX::Abstract::Vector&
	update (	double alpha,
				const NOX::Abstract::Vector& a,
				double beta,
				const NOX::Abstract::Vector& b,
				double gamma=0.0)
	{ return update(alpha,dynamic_cast<const NOXVector<T_DOF>&>(a),beta,dynamic_cast<const NOXVector<T_DOF>&>(b),gamma); }

	static const unsigned int
	dof ()
	{ return DOF; }

	NOXVector<DOF>
	toNOXVector ( ) const
	{ return *this; }
};

#endif
