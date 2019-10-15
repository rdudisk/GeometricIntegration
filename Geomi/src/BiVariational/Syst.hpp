#ifndef DEF_BIVARIATIONAL_SYST
#define DEF_BIVARIATIONAL_SYST

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

namespace BiVariational {

template <typename T_SCALAR, typename T_Q, typename T_VEL>
class Syst
{
	typedef		std::vector<T_SCALAR>		CoordVec;
	typedef		std::vector<T_VEL>			VelVec;
protected:
	CoordVec	m_coord;
	T_Q			m_pos;
	VelVec		m_vel;

public:
	Syst ( )
	{ }

	Syst (CoordVec coord)
	: m_coord(coord)
	{ }

	Syst (CoordVec coord, T_Q pos)
	: m_coord(coord), m_pos(pos)
	{ }

	Syst (CoordVec coord, T_Q pos, VelVec vel)
	: m_coord(coord), m_pos(pos), m_vel(vel)
	{ }

	~Syst ( )
	{ }

	CoordVec
	coord (void) const
	{ return m_coord; }

	void
	coord (CoordVec b)
	{ m_coord = b; }

	void
	coord (int i, T_SCALAR s)
	{ m_coord[i] = s; }
	
	T_Q
	pos ( ) const
	{ return m_pos; }
	
	void
	pos (T_Q p)
	{ m_pos = p; }

	VelVec
	vel ( ) const
	{ return m_vel; }

	void
	vel (VelVec v)
	{ m_vel = v; }
};

} // namespace BiVariational

#endif
