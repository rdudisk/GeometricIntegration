#ifndef DEF_BIVARIATIONAL_SOLVEME
#define DEF_BIVARIATIONAL_SOLVEME

namespace BiVariational {

template <typename T_SCALAR, typename T_Q, typename T_VEL>
class StepInternals : public NOX::Epetra::Interface::Required
{
private:
	Abstract::LieProblem<T_SCALAR,T_Q,T_VEL>* m_problem;
	Epetra_Map* m_map;
	Epetra_Comm* m_comm;

	Teuchos::RCP<Epetra_Vector> m_initialGuess;

	int m_numLocalSystems;
	int m_time_index;

	std::vector<Eigen::Matrix<T_SCALAR,6,1>> lambda;
	std::vector<Eigen::Matrix<T_SCALAR,6,1>> mu;

public:
	StepInternals (Abstract::LieProblem<T_SCALAR,T_Q,T_VEL>& problem, Epetra_Comm& comm)
	:	m_problem(&problem),
		m_comm(&comm)
	{
		int i,j;
		int dof = 6; //m_problem->DOF;
		
		int myPID = m_comm->MyPID();
		int numProcs = m_comm->NumProc();

		size_t n_space_steps = m_problem->size(1);
		const int numGlobalEntries = dof * n_space_steps;

		const int indexBase = 0;

		m_numLocalSystems = (n_space_steps+numProcs-myPID-1)/numProcs;
		int numLocalEntries = m_numLocalSystems*dof;
		int* gblIndList = new int [numLocalEntries];

		for (i=0; i<m_numLocalSystems; i++) {
			for (j=0; j<dof; j++) {
				gblIndList[i*dof+j] = myPID*dof + i*dof*numProcs + j;
			}
		}

		m_map = new Epetra_Map(numGlobalEntries,numLocalEntries,gblIndList,indexBase,comm);

		if (gblIndList != NULL) {
			delete [] gblIndList;
			gblIndList = NULL;
		}

		m_initialGuess = Teuchos::rcp(new Epetra_Vector(*m_map));

		m_time_index = 0;

		lambda.resize(n_space_steps);
		mu.resize(n_space_steps);
	}

	~StepInternals ()
	{
		delete m_map;
	}

	Teuchos::RCP<Epetra_Vector>
	getInitialGuess ()
	{
		int i,j;
		int* gblIndList = m_map->MyGlobalElements();
		int numLocalEntries = m_map->NumMyElements();
		int space_index, time_index;
		int dof = 6;
		if (numLocalEntries != dof*m_numLocalSystems)
			std::cout << "Error: numLocalEntries != dof*m_numLocalSystems, " << numLocalEntries << " != " << dof << "*" << m_numLocalSystems << std::endl;

		T_VEL vel;

		for (j=0; j<m_numLocalSystems; j++) {
			space_index = gblIndList[j*dof]/dof;
			if (m_time_index>0)
				time_index = m_time_index-1;
			else
				time_index = 0;
			vel = m_problem->vel(time_index,space_index)[0];

			for (i=0; i<dof; i++) {
				(*m_initialGuess)[j*dof+i] = vel[i];
			}
		}

		return m_initialGuess;
	}

	bool
	initLimit ( )
	{
		// On suppose que la configuration du systeme est donné pour t = 0,h
		int j;
		double h;
		T_VEL xi;

		for(j=0; j<m_problem->size(1); j++) {
			h = m_problem->step_size(0,j,0);
			xi = (1.0/h)*T_VEL::cay_inv(m_problem->pos(0,j).inverse()*m_problem->pos(1,j));
			mu[j] = (h*xi).dCayRInv().transpose()*m_problem->dLdv1(xi);
			std::cout << "h: " << h << std::endl;
			std::cout << "fkk+1: " << std::endl << (m_problem->pos(0,j).inverse()*m_problem->pos(1,j)).matrix() << std::endl;
			std::cout << "xi:" << std::endl << xi << std::endl;
			std::cout << "mu [" << j << "]: " << std::endl << mu[j] << std::endl;
		}

		return true;
	}

	bool
	setupStep ( )
	{
		// TODO: parallel
		
		int i = m_time_index, j;
		double h = m_problem->step_size(i,0,0);
		double s;

		Eigen::Matrix<double,6,1> mu_prev;
		T_VEL eta, xi_prev;

		for (j=0; j<m_problem->size(1)-1; j++) {
			s = m_problem->step_size(i,j,1);

			eta = ((1.0/s)*T_VEL::cay_inv(m_problem->pos(i,j).inverse()*(m_problem->pos(i,j+1)))).toVector();
			xi_prev = (1.0/h)*T_VEL::cay_inv(m_problem->pos(i-1,j).inverse()*m_problem->pos(i,j));
			lambda[j] = (s*eta).dCayRInv().transpose()*m_problem->dLdv1(eta);
			mu_prev = mu[j];
			mu[j] = T_VEL(mu_prev).Ad_star(m_problem->pos(i-1,j).inverse()*m_problem->pos(i,j)).toVector() + (h/s)*lambda[j];
			if (j>0)
				mu[j] += (-(h/s)*T_VEL(lambda[j-1]).Ad_star(m_problem->pos(i,j-1).inverse()*m_problem->pos(i,j))).toVector();
		}

		return true;
	}

	bool
	computeF (const Epetra_Vector& x, Epetra_Vector& f, const FillType flag = Residual)
	{
		int i,j;
		int* gblIndList = m_map->MyGlobalElements();

		const int dof = 6;

		Eigen::Matrix<double,dof,1> tmp_f;
		Eigen::Matrix<double,dof,1> tmp_x;
		T_VEL xi;

		int time_index = m_time_index, space_index;
		double s, h = m_problem->step_size(time_index,0,0);

		for (j=0; j<m_numLocalSystems; j++) {
			space_index = gblIndList[j*dof]/dof;

			for (i=0; i<dof; i++)
				tmp_x[i] = x[i+j*dof];
			xi = T_VEL(tmp_x);

			tmp_f =	(h*xi).dCayRInv().transpose()*m_problem->dLdv0(xi) - mu[space_index];

			for (i=0; i<dof; i++) {
				f[j*dof+i] = tmp_f[i];
			}
			//std::cout << tmp_f << std::endl;
		}

		return true;
	}

	/*
	void
	test ()
	{
		int i;
		int* gblIndList = m_map->MyGlobalElements();
		int numLocalEntries = m_map->NumMyElements();
		int myPID = m_comm->MyPID();

		Epetra_Vector f(*m_map);
		
		for (i=0; i<numLocalEntries; i++) {
			f[i] = 1.0*(1+myPID);
		}

		std::cout << f << std::endl;
	}
	*/
};

} // namespace BiVariational

#endif
