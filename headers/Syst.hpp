#ifndef DEF_SYST
#define DEF_SYST

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "csv.hpp"

template <typename M, typename Q>
class Syst
{
protected:
	M m_base; /**< Base space variable */
	Q m_pos; /**< Configuration space variable */

public:
	Syst<M,Q> ( )
	{ }

	Syst<M,Q> (M base)
	: m_base(base)
	{ }

	Syst<M,Q> (M base, Q posl)
	: m_base(base), m_pos(pos)
	{ }

	~Syst<M,Q> ( )
	{ }

	M
	base (void) const
	{
		return m_base;
	}

	void
	base(M _base)
	{
		m_base = _base;
	}
	
	Q
	pos ( ) const
	{
		return m_pos;
	}
	
	void
	pos (Q _pos)
	{
		m_pos = _pos;
	}

	friend std::ostream&
	operator<< (std::ostream& stream, const Syst<M,Q>& s)
	{
		stream << "(" << s.m_base << "," << s.m_pos << ")";
		return stream;
	}

	std::string
	csvString (const std::string sep)
	{
		std::ostringstream ss;
		ss << ::csvString<M>(m_base,sep) << sep << ::csvString<Q>(m_pos,sep);
		std::string str(ss.str());
		return str;
	}
};

#endif
