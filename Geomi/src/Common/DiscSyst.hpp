#ifndef DEF_COMMON_DISCSYST
#define DEF_COMMON_DISCSYST

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

//#include "Common_Syst.hpp"

/**
 * Class for discrete jet space.
 * A discrete jet space consists of an array over the base space of jet space elements.
 * \tparam M the class representing the base space.
 * \tparam Q the class representing the configuration space.
 */
template <typename T_M, typename T_Q>
class DiscSyst
{
protected:
	std::vector<Syst<T_M,T_Q>> m_node;

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
	~DiscSyst<T_M,T_Q> ( )
	{ }

	Syst<T_M,T_Q> const&
	operator[] (size_t index) const
	{
		if (index<0 || index>=m_node.size()) {
			throw std::out_of_range("Error : DiscSyst[] index out of range");
		}
		return m_node[index];
	}

	Syst<T_M,T_Q>&
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

	/* Accessors and setters */

	T_M
	base (const size_t& i) const
	{
		return m_node[i].base();
	}

	T_Q
	pos (const size_t& i) const
	{
		return m_node[i].pos();
	}

	void
	base (const size_t& i, const T_M& b)
	{
		m_node[i].base(b);
	}

	void
	pos (const size_t& i, const T_Q& p)
	{
		m_node[i].pos(p);
	}

	static const unsigned int
	dof ( )
	{
		return T_Q::DOF;
	}

	/**
	 * Set the base space discretization as a linear interpolation consisting in \p n_steps steps
	 * between \p inf_lim and \p sup_lim.
	 * \warning This completely overwrites the entire data content of `*this`.
	 */
	void
	baselinspace (const T_M& inf_lim, const T_M& sup_lim, const size_t& n_steps)
	{
		if (n_steps<0)
			throw std::domain_error("Error : DiscSyst.linspace size argument must be null or positive");
		m_node = std::vector<Syst<T_M,T_Q>>();
		if (n_steps>0) {
			float fn_steps = (float)n_steps;
			for (float i=0.0f; i<fn_steps; i+=1.0f) {
				m_node.push_back(Syst<T_M,T_Q>(inf_lim + (sup_lim-inf_lim)*i/(fn_steps-1.0f)));
			}
		}
	}

	/**
	 * Set the base space discretization as a linear interpolation constisting in \p n_steps steps of length
	 * \p step_size starting at \p inf_lim.
	 * \warning This completely overwrites the entire data content of `*this`.
	 */
	void
	baselinstep (T_M inf_lim, T_M step_size, size_t n_steps)
	{
		if (n_steps<0)
			throw std::domain_error("Error : DiscSyst.linstep size argument must be null or positive");
		m_node = std::vector<Syst<T_M,T_Q>>();
		if (n_steps>0) {
			for (int i=0; i<n_steps; i++) {
				m_node.push_back(Syst<T_M,T_Q>(inf_lim+i*step_size));
			}
		}
	}

	friend std::ostream&
	operator<< (std::ostream& stream, const DiscSyst<T_M,T_Q>& ds)
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
