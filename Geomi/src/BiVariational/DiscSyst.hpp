#ifndef DEF_BIVARIATIONAL_DISCSYST
#define DEF_BIVARIATIONAL_DISCSYST

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

template <typename T_SCALAR, typename T_Q, typename T_VEL>
class DiscSyst
{
	typedef		std::vector<T_SCALAR>		CoordVec;
	typedef		std::vector<T_VEL>			VelVec;
protected:
	std::vector<Syst<T_SCALAR,T_Q,T_VEL>> m_node;
	size_t m_size[2];

public:
	DiscSyst ( )
	{
		m_size[0] = 0;
		m_size[1] = 0;
		m_node.clear();
	}

	~DiscSyst ( )
	{ }

	const size_t
	getIndex (const size_t i, const size_t j) const
	// TODO: check range return value
	{ return i*m_size[1]+j; }

	/* Accessors and setters */

	Syst<T_SCALAR,T_Q,T_VEL>
	node (const size_t i, const size_t j) const
	{ return m_node[getIndex(i,j)]; }

	size_t
	size (size_t i) const
	// TODO: check i range
	{ return m_size[i]; }

	CoordVec
	coord (const size_t i, const size_t j) const
	{ return m_node[getIndex(i,j)].coord(); }

	void
	coord (const size_t i, const size_t j, const CoordVec b)
	{ m_node[getIndex(i,j)].coord(b); }

	T_Q
	pos (const size_t i) const
	{ return m_node[getIndex(i,j)].pos(); }

	void
	pos (const size_t i, const size_t j, const T_Q p)
	{ m_node[getIndex(i,j)].pos(p); }

	VelVec
	vel (const size_t i, const size_t j) const
	{ return m_node[getIndex(i,j)].vel(); }

	void
	vel (const size_t i, const size_t j, const VelVec v)
	{ m_node[getIndex(i,j)].vel(v); }

	static const unsigned int
	dof ( )
	{ return T_Q::DOF; }
};

#endif
