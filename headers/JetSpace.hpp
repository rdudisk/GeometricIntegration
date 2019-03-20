#ifndef DEF_JET_SPACE
#define DEF_JET_SPACE

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

#include "csv.hpp"

/** Petite description de test.
 *  Sur une seconde ligne.
 *  Bla.
 */
template <typename M, typename Q, typename TQ>
class JetSpace
{
protected:
	M m_base; /**< Base space variable */
	Q m_pos; /**< Configuration space variable */
	TQ m_vel; /**< Velocity space variable */

public:
	JetSpace<M,Q,TQ> ( )
	{ }

	JetSpace<M,Q,TQ> (M base)
	: m_base(base)
	{ }

	JetSpace<M,Q,TQ> (M base, Q pos, TQ vel)
	: m_base(base), m_pos(pos), m_vel(vel)
	{ }

	~JetSpace<M,Q,TQ> ( )
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

	TQ
	vel ( ) const
	{
		return m_vel;
	}

	void
	vel (TQ _vel)
	{
		m_vel = _vel;
	}

	friend std::ostream&
	operator<< (std::ostream& stream, const JetSpace<M,Q,TQ>& js)
	{
		stream << "(" << js.m_base << "," << js.m_pos << "," << js.m_vel << ")";
		return stream;
	}

	std::string
	csvString (const std::string sep)
	{
		std::ostringstream ss;
		/*if (std::is_fundamental<M>::value) {
			ss << m_base;
		} else {
			ss << m_base.csvString(sep);
		}
		ss << sep;
		if (std::is_fundamental<Q>::value) {
			ss << m_pos;
		} else {
			ss << m_pos.csvString(sep);
		}
		ss << sep;
		if (std::is_fundamental<TQ>::value) {
			ss << m_vel;
		} else {
			ss << m_vel.csvString(sep);
		}
		*/
		//ss << m_base << sep << m_pos << sep << m_vel;
		ss << ::csvString<M>(m_base,sep) << sep << ::csvString<Q>(m_pos,sep) << sep << ::csvString<TQ>(m_vel,sep);
		std::string str(ss.str());
		return str;
		//return m_base.toString() + separator + m_pos.toString() + separator + m_vel.toString();
	}
};

#endif
