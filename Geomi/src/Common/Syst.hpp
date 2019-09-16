#ifndef DEF_COMMON_SYST
#define DEF_COMMON_SYST

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

//#include "Common_csv.hpp"

template <typename T_M, typename T_Q>
class Syst
{
protected:
	T_M m_base; /**< Base space variable */
	T_Q m_pos; /**< Configuration space variable */

public:
	Syst<T_M,T_Q> ( )
	{ }

	Syst<T_M,T_Q> (T_M base)
	: m_base(base)
	{ }

	Syst<T_M,T_Q> (T_M base, T_Q pos)
	: m_base(base), m_pos(pos)
	{ }

	~Syst<T_M,T_Q> ( )
	{ }

	T_M
	base (void) const
	{ return m_base; }

	void
	base(T_M b)
	{ m_base = b; }
	
	T_Q
	pos ( ) const
	{ return m_pos; }
	
	void
	pos (T_Q p)
	{ m_pos = p; }

	friend std::ostream&
	operator<< (std::ostream& stream, const Syst<T_M,T_Q>& s)
	{
		stream << "(" << s.m_base << "," << s.m_pos << ")";
		return stream;
	}

	std::string
	csvString (const std::string sep)
	{
		std::ostringstream ss;
		ss << ::csvString<T_M>(m_base,sep) << sep << ::csvString<T_Q>(m_pos,sep);
		std::string str(ss.str());
		return str;
	}
};

#endif
