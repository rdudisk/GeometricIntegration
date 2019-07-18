#ifndef DEF_PROBLEM
#define DEF_PROBLEM

/*
template <typename M, typename Q, typename Alg>
class Problem
{
private:
	// The actual physical configurations
	DiscSyst<M,Q>& m_ds;

	ProblemInterface& m_interface;

public:

	Problem (DiscSyst<M,Q>& ds, ProblemInterface& interface) :
		m_ds(ds),
		m_interface(interface)
	{ }

	~Problem ()
	{ }

	bool
	computeA (Alg& A, const Q& x)
	{
		return m_interface.computeA(A,y);
	}

	bool
	computeJacobianA (std::vector<Alg>& JA, const Q& x)
	{
		return m_interface.computeJacobianA(JA,y);
	}

};
*/

#endif

