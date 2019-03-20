#ifndef DEF_VEC
#define DEF_VEC

#include <Eigen/Dense>

#include <gsl/gsl_vector.h>

template <typename T,size_t N>
class Vec : public Eigen::Matrix<T,N,1> {
public:
	//using Eigen::Matrix<T,N,1>::Matrix;
	
	Vec<T,N> ( )
	: Eigen::Matrix<T,N,1>()
	{ }

	/* see Eigen : inheriting from matrix for the two following methods */
	template<typename OtherDerived>
	Vec<T,N> (const Eigen::MatrixBase<OtherDerived>& other)
	: Eigen::Matrix<T,N,1>(other)
	{ }

	template<typename OtherDerived>
	Vec<T,N>&
	operator= (const Eigen::MatrixBase<OtherDerived>& other)
	{
		this->Eigen::Matrix<T,N,1>::operator=(other);
		return *this;
	}

	static size_t
	dof ( )
	{
		return N;
	}

	static Vec<T,N>
	cast_from_gsl_vector (const gsl_vector* x)
	{
		Vec<T,N> res;
		for (int i=0; i<N; i++) {
			res(i) = gsl_vector_get(x,i);
		}
		return res;
	}

	gsl_vector*
	cast_gsl_vector ( )
	{
		gsl_vector *x = gsl_vector_alloc(N);
		for (int i=0; i<N; i++) {
			gsl_vector_set(x,i,(double)(*this)(i));
		}
		return x;
	}

	// Operateurs

	/*
	using Eigen::Matrix<T,N,1>::operator+=;

	Vec<T,N> operator+ (const Vec<T,N>& v) const {
		Vec<T,N> res(*this);
		res += v;
		return res;
	}
	
	using Eigen::Matrix<T,N,1>::operator*=;
	using Eigen::Matrix<T,N,1>::operator*;

	using Eigen::Matrix<T,N,1>::operator/=;
	using Eigen::Matrix<T,N,1>::operator/;
	*/
};

#endif
