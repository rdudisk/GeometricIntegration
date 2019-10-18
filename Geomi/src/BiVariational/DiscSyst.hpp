#ifndef DEF_BIVARIATIONAL_DISCSYST
#define DEF_BIVARIATIONAL_DISCSYST

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

namespace BiVariational {

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

	void
	setSize (const size_t i, const size_t j)
	{
		m_size[0] = i;
		m_size[1] = j;
		m_node.resize(i*j);

		std::vector<T_SCALAR> default_coord(2);
		default_coord[0] = 0.0; default_coord[1] = 0.0;
		T_Q default_pos;
		std::vector<T_VEL> default_vel(2);
		default_vel[0] = T_VEL(); default_vel[1] = T_VEL();

		Syst<T_SCALAR,T_Q,T_VEL> default_node(default_coord,default_pos,default_vel);

		int k;

		for (k=0; k<i*j; k++) {
			m_node[k] = default_node;
		}
	}

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

	T_SCALAR
	step_size (const size_t i, const size_t j, const size_t dir)
	{
		T_SCALAR t0 = m_node[getIndex(i,j)].coord()[dir];
		T_SCALAR t1;
		if ((dir==0 && i==m_size[0]-1) || (dir==1 && j==m_size[1]-1)) {
			std::cout << "Error, I can't give you the step in dir " << dir << " at point " << i << "," << j << std::endl;
			return 0;
		}
		if (dir==0)
			t1 = m_node[getIndex(i+1,j)].coord()[0];
		else
			t1 = m_node[getIndex(i,j+1)].coord()[1];
		return t1-t0;
	}

	T_Q
	pos (const size_t i, const size_t j) const
	{ return m_node[getIndex(i,j)].pos(); }

	void
	pos (const size_t i, const size_t j, const T_Q p)
	{ m_node[getIndex(i,j)].pos(p); }

	VelVec
	vel (const size_t i, const size_t j) const
	{ return m_node[getIndex(i,j)].vel(); }

	T_VEL
	vel (const size_t i, const size_t j, const size_t dir) const
	{ return m_node[getIndex(i,j)].vel(dir); }

	void
	vel (const size_t i, const size_t j, const VelVec v)
	{ m_node[getIndex(i,j)].vel(v); }

	void
	vel (const size_t i, const size_t j, const size_t dir, const T_VEL v)
	{ m_node[getIndex(i,j)].vel(dir,v); }

	static const unsigned int
	dof ( )
	{ return T_Q::DOF; }

	/**
	 * Set the base space discretization as a linear interpolation constisting in \p n_steps steps of length
	 * \p step_size starting at \p inf_lim.
	 * \warning This completely overwrites the entire data content of `*this`.
	 */
	void
	baselinstep (T_SCALAR time_inf_lim, T_SCALAR time_step_size, T_SCALAR space_inf_lim, T_SCALAR space_step_size)
	{
		int i,j;
		std::vector<T_SCALAR> c(2);
		for (i=0; i<m_size[0]; i++) {
			c[0] = time_inf_lim+i*time_step_size;
			for (j=0; j<m_size[1]; j++) {
				c[1] = space_inf_lim+j*space_step_size;
				coord(i,j,c);
			}
		}
	}

};

}

#endif
