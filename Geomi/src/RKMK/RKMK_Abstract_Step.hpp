#ifndef DEF_RKMK_ABSTRACT_STEP
#define DEF_RKMK_ABSTRACT_STEP

#include <algorithm>	// std::max
#include <iterator>		// C style array copy


namespace RKMK {
namespace Abstract {
/**
 * Abstract template class for a RKMK step.
 * This class is used by RKMK::ExplicitStep and RKMK::DiagonalStep.
 */

template <typename T_M,
		  typename T_Q,
		  typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES>
class Step
{ 
protected:
	RKMK::StepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>* m_internals;
	RKMK::Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& m_interface;
	char m_type;

public:
	static const char TYPE_UNKNOWN = 0;
	static const char TYPE_EXPLICIT = 1;
	static const char TYPE_DIAGONAL_IMPLICIT = 2;

public:
	Step<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (RKMK::Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface)
	//													 Abstract::RKMKStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>& internals)
	:	m_interface(interface)
		//m_internals(internals)
		//m_internals(Abstract::RKMKStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>::newFromProblem(interface))
	{
		m_internals = new RKMK::StepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>(interface);
		setType();
	}

	/*
	RKMKStep<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> (Abstract::Problem<T_M,T_Q,T_LIE_ALGEBRA>& interface,
														 std::vector<double> a,
														 std::vector<double> b)
	:	m_interface(interface),
		m_internals(Abstract::RKMKStepInternals<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES>::newFromProblem(interface))
	{ m_internals.setCoeffs(a,b); }
	*/

	~Step<T_M,T_Q,T_LIE_ALGEBRA,T_N_INTERNAL_STAGES> ()
	{ }

	bool
	setCoeffs (std::vector<double> a, std::vector<double> b) const
	{ return m_internals->setCoeffs(a,b); }

	void
	setData (T_M h_var, T_Q y0_var)
	{ m_internals->setData(h_var,y0_var); }

	void
	setType ( )
	{
		int i,j;
		bool tmp_isZero;
		bool isExplicit = true;
		bool isDIRK = true;

		for (i=0; i<T_N_INTERNAL_STAGES; i++) {
			isExplicit &= isZero<double>(m_internals->a_coeffs(i,i));
			for (j=i+1; j<T_N_INTERNAL_STAGES; j++) {
				tmp_isZero = isZero<double>(m_internals->a_coeffs(i,j));
				isExplicit &= tmp_isZero;
				isDIRK &= tmp_isZero;
			}
		}	

		if (isExplicit)
			m_type = TYPE_EXPLICIT;
		else if (isDIRK)
			m_type = TYPE_DIAGONAL_IMPLICIT;
		else
			m_type = TYPE_UNKNOWN;
	}

	const char
	type ()
	{ return m_type; }

	virtual const T_Q
	makeStep () = 0;
};
} // namespace Abstract
} // namespace RKMK

#endif
