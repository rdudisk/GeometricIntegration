#ifndef DEF_DISCRETE_SYSTEM_SPACE
#define DEF_DISCRETE_SYSTEM_SPACE

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "Syst.hpp"

/**
 * Class for discrete jet space.
 * A discrete jet space consists of an array over the base space of jet space elements.
 * \tparam M the class representing the base space.
 * \tparam Q the class representing the configuration space.
 */
template <typename M, typename Q>
class DiscSyst
{
protected:
	std::vector<Syst<M,Q>> m_node;

	int
	__write2csv__ (const std::string filename, std::ios_base::openmode mode, const std::string sep)
	{
		std::ofstream file;
		file.open(filename,mode);
		if (file.is_open()) {
			for (int i=0; i<size(); i++) {
				file << m_node[i].csvString(sep) << std::endl;
			}
			file.close();
			return 0;
		} else {
			return 1;
		}
	}

public:
	DiscSyst<M,Q> ( )
	{ }

	~DiscSyst<M,Q> ( )
	{ }

	Syst<M,Q> const&
	operator[] (size_t index) const
	{
		if (index<0 || index>=m_node.size()) {
			throw std::out_of_range("Error : DiscSyst[] index out of range");
		}
		return m_node[index];
	}

	Syst<M,Q>&
	operator[] (size_t index)
	{
		if (index<0 || index>=m_node.size()) {
			throw std::out_of_range("Error : DiscSyst[] index out of range");
		}
		return m_node[index];
	}

	size_t
	size() const
	{
		return m_node.size();
	}

	/* Accdessors and setters */

	M
	base (const size_t& i) const
	{
		return m_node[i].base();
	}

	Q
	pos (const size_t& i) const
	{
		return m_node[i].pos();
	}

	void
	base (const size_t& i, const M& _base)
	{
		m_node[i].base(_base);
	}

	void
	pos (const size_t& i, const Q& _pos)
	{
		m_node[i].pos(_pos);
	}

	size_t
	dof ( )
	{
		return Q::dof();
	}

	/**
	 * Set the base space discretization as a linear interpolation consisting in \p n_steps steps
	 * between \p inf_lim and \p sup_lim.
	 * \warning This completely overwrites the entire data content of `*this`.
	 */
	void
	baselinspace (const M& inf_lim, const M& sup_lim, const size_t& n_steps)
	{
		if (n_steps<0)
			throw std::domain_error("Error : DiscSyst.linspace size argument must be null or positive");
		m_node = std::vector<Syst<M,Q>>();
		if (n_steps>0) {
			float fn_steps = (float)n_steps;
			for (float i=0.0f; i<fn_steps; i+=1.0f) {
				m_node.push_back(Syst<M,Q>(inf_lim + (sup_lim-inf_lim)*i/(fn_steps-1.0f)));
			}
		}
	}

	/**
	 * Set the base space discretization as a linear interpolation constisting in \p n_steps steps of length
	 * \p step_size starting at \p inf_lim.
	 * \warning This completely overwrites the entire data content of `*this`.
	 */
	void
	baselinstep (M inf_lim, M step_size, size_t n_steps)
	{
		if (n_steps<0)
			throw std::domain_error("Error : DiscSyst.linstep size argument must be null or positive");
		m_node = std::vector<Syst<M,Q>>();
		if (n_steps>0) {
			for (int i=0; i<n_steps; i++) {
				m_node.push_back(Syst<M,Q>(inf_lim+i*step_size));
			}
		}
	}

	friend std::ostream&
	operator<< (std::ostream& stream, const DiscSyst<M,Q>& ds)
	{
		for (int i=0; i<ds.size(); i++) {
			stream << ds[i];
			if (i!=ds.size())
				stream << "\t";
			else
				stream << std::endl;
		}
		return stream;
	}

	int
	append2csv (const std::string filename, const std::string sep=",")
	{
		return __write2csv__(filename,std::ios_base::app,sep);
	}

	int
	write2csv (const std::string filename, const std::string sep=",")
	{
		return __write2csv__(filename,std::ios_base::trunc,sep);
	}
};

#endif
