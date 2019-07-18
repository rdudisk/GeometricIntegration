#ifndef NOX_MULTIVECTOR
#define NOX_MULTIVECTOR

/*
#include "NOX_Common.H"
#include "NOX_Abstract_MultiVector.H"

#include <Eigen/Base>

#include "NOXVector.hpp"

template <typename T, unsigned int dof>
class NOXMultiVector : public virtual NOX::Abstract::MultiVector
{
public:
	~NOXMultiVector<T,dof> ()
	{ }

	NOXVector<T,dof> ()
	{ }

	template <typename OtherDerived>
	NOXVector<T,dof> (const Eigen::MatrixBase<OtherDerived>& other) :
		Eigen::Matrix<T,dof,1>(other)
	{ }

	template <typename OtherDerived>
	NOX::Abstract::Vector&
	operator= (const Eigen::MatrixBase<OtherDerived>& other)
	{
		this->Eigen::Matrix<T,dof,1>::operator=(other);
		return *this;
	}

	// Copy constructor with NOX::CopyType : just ignore type
	NOXVector<T,dof> (const NOXVector& source, NOX::CopyType type) :
		NOXVector<T,dof>(source)
	{ }

	NOX::Abstract::Vector&
	abs (const NOXVector& y)
	{
		this = y.Eigen::Matrix<T,dof>::abs();
		return *this;
	}

	NOX::Abstract::Vector&
	abs (const NOX::Abstract::Vector& y)
	{
		return abs(dynamic_cast<const NOXVector<T,dof>&>(y));
	}

	Teuchos::RCP<NOX::Abstract::Vector>
	clone (NOX::CopyType type=NOX::DeepCopy) const
	{
		Teuchos::RCP<NOX::Abstract::Vector> tmp;
		tmp = Teuchos::rcp(new NOXVector(*this,type));
		return tmp;
	}

	// TODO:
	// override createMultiVector()
	// pour ca il faut une implementation de NOX::Abstract::MultiVector, qui par defaut (le cas ici) est NOX::MultiVector

	double
	innerProduct (const NOXVector<T,dof>& y) const
	{
		return (this.dot(y).sum());
	}

	double
	innerProduct (const NOX::Abstract::Vector& y) const
	{
		return innerProduct(dynamic_cast<const NOXVector<T,dof>&>(y));
	}

	NOX::Abstract::Vector&
	init (double gamma)
	{
		this.setConstant(gamma);
		return *this;
	}

	NOX::size_type
	length () const
	{
		return dof;
	}

	double
	norm (NOX::Abstract::Vector::NormType type=NOX::Abstract::Vector::TwoNorm) const
	{
		double result;

		switch (type) {
			case NOX::Abstract::Vector::MaxNorm:
				result = this.maxCoeff();
				break;
			case NOX::Abstract::Vector::OneNorm:
				result = this.lpNorm<1>();
			case NOX::Abstract::Vector::TwoNorm:
				result = this.norm();
				break;
			default:
				break;
		}
		
		return result;
	}

	double
	norm (const NOXVector& weights) const
	{
		return sqrt(((this.cwiseProduct(this)).cwiseProduct(weights)).sum());
	}

	double
	norm (const NOX::Abstract::Vector& weights) const
	{
		return norm(dynamic_cast<const NOXVector<T,dof>&>(weights));
	}

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
			NOX::Random:setSeed(seed);
		
		for (int i=0; i<dof; i++)
			this[i] = NOX::Random::number();

		return *this;
	}

	// TODO: definir operator=(const NOXVector&)
	NOX::Abstract::Vector&
	operator= (const NOX::Abstract::Vector& y)
	{
		return operator=(dynamic_cast<const NOXVector<T,dof>&>(y));
	}

	NOX::Abstract::Vector&
	reciprocal (const NOXVector& y)
	{
		this = this.cwiseInverse();
		return *this;
	}

	NOX::Abstract::Vector&
	reciprocal (const NOX::Abstract::Vector& y)
	{
		return reciprocal(dynamic_cast<const NOXVector<T,dof>&>(y));
	}

	NOX::Abstract::Vector&
	scale (double gamma)
	{
		this = this*gamma;
		return *this;
	}

	NOX::Abstract::Vector&
	scale (const NOXVector<T,dof>& a)
	{
		this = this.cwiseProduct(a);
		return *this;
	}

	NOX::Abstract::Vector&
	scale (const NOX::Abstract::Vector& a)
	{
		return scale(dynamic_cast<const NOXVector<T,dof>&>(a));
	}

	NOX::Abstract::Vector&
	update (double alpha, const NOXVector<T,dof>& a, double gamma=0.0)
	{
		this = alpha*a + this*gamma;
		return *this;
	}

	NOX::Abstract::Vector&
	update (double alpha, const NOX::Abstract::Vector& a, double gamma=0.0)
	{
		return update(alpha,dynamic_cast<const NOXVector<T,dof>&>(a),gamma);
	}

	NOX::Abstract::Vector&
	update (	double alpha,
				const NOXVector<T,dof>& a,
				double beta,
				const NOXVector<T,dof>& b,
				double gamma=0.0)
	{
		this = alpha*a + beta*b + this*gamma;
		return *this;
	}

	NOX::Abstract::Vector&
	update (	double alpha,
				const NOX::Abstract::Vector& a,
				double beta,
				const NOX::Abstract::Vector& b,
				double gamma=0.0)
	{
		return update(alpha,dynamic_cast<const NOXVector<T,dof>&>(a),beta,dynamic_cast<const NOXVector<T,dof>&>(b),gamma);
	}
};
*/

#endif

