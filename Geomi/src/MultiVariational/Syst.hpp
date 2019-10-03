#ifndef DEF_MULTIVARIATIONAL_SYST
#define DEF_MULTIVARIATIONAL_SYST

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

namespace MultiVariational {

template <typename T_Q, int T_N_INDEP>
class Syst
{
protected:
	std::vector<double> m_base;
	T_Q m_pos;

public:
	Syst<T_Q,T_N_INDEP> ( )
	{ }

	Syst<T_Q,T_N_INDEP> (std::vector<double> base)
	: m_base(base)
	{ }

	Syst<T_Q,T_N_INDEP> (std::vector<double> base, T_Q pos)
	: m_base(base), m_pos(pos)
	{ }

	~Syst<T_Q,T_N_INDEP> ( )
	{ }

	std::vector<double>
	base (void) const
	{ return m_base; }

	void
	base(std::vector<double> b)
	{ m_base = b; }
	
	T_Q
	pos ( ) const
	{ return m_pos; }
	
	void
	pos (T_Q p)
	{ m_pos = p; }

	std::string
	csvString (const std::string sep)
	{
		std::ostringstream ss;
		ss << ::csvString<std::vector<double>>(m_base,sep) << sep << ::csvString<T_Q>(m_pos,sep);
		std::string str(ss.str());
		return str;
	}
};

} // namespace MultiVariational

#endif
