#ifndef DEF_COMMON_DISCSYST
#define DEF_COMMON_DISCSYST

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

namespace MultiVariational {

template <typename T_Q, int T_N_INDEP>
class DiscSyst
{
protected:
	std::vector<MltiVariational::Syst<T_Q,T_N_INDEP>> m_node;
	std::vector<int> m_size; /** size along each dimension */
	int m_total_size;

	int
	__write2csv__ (const std::string filename, std::ios_base::openmode mode, const std::string sep)
	{
		// TODO
		return 0;
		/*std::ofstream file;
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
		*/
	}

public:
	~DiscSyst<T_Q,T_N_INDEP> (std::vector<int> size)
	{
		// TODO: check size
		m_size = size;
		m_totat_size = 1;
		for (int i=0; i<T_N_INDEP; i++)
			m_total_size *= m_size[i];
		m_node.resize(m_total_size);
	}

	std::vector<int>
	size() const
	{ return m_size; }

	/* Accessors and setters */

	int
	getIndex (const int[] index) const
	{
		int ret = 0;
		for (int i=0; i<T_N_INDEP; i++)
			ret = ret*m_size[i] + index[i];
		return ret;
	}

	std::vector<double>
	base (const int[] index) const
	{ return m_node[getIndex(i)].base(); }

	T_Q
	pos (const int[] index) const
	{ return m_node[getIndex(i)].pos(); }

	void
	base (const int[] index, const std::vector<double>& b)
	{ m_node[getIndex(i)].base(b); }

	void
	pos (const int[] index, const T_Q& p)
	{ m_node[getIndex(i)].pos(p); }

	static const unsigned int
	dof ( )
	{ return T_Q::DOF; }

	/**
	 * Set the base space discretization as a linear interpolation constisting in \p n_steps steps of length
	 * \p step_size starting at \p inf_lim.
	 * \warning This completely overwrites the entire data content of `*this`.
	 */
	void
	baselinstep (std::vector<double> init_date,std::vector<double> step_size, std::vector<int> n_steps);

	int
	append2csv (const std::string filename, const std::string sep=",")
	{ return __write2csv__(filename,std::ios_base::app,sep); }

	int
	write2csv (const std::string filename, const std::string sep=",")
	{ return __write2csv__(filename,std::ios_base::trunc,sep); }
};

} // namespace MultiVariational

template <T_Q>
MultiVariational::DiscSyst<T_Q,2>::baselinstep (std::vector<double> init_date, std::vector<double> step_size, std::vector<int> n_steps)
{
	std::vector<double> date = init_date;
	int i,j;
	for (i=0; n<n_steps[0]; i++) {
		date[0] = init_date
		for (j=0; j<n_steps[1]; j++) {
			date[1] = init_date[1]+j*step_size[1];
			m_node.base(getIndex(i,j),Syst<T_Q,2>(date));
		}
	}
}

#endif
