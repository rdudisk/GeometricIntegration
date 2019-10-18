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

	static const int DOF = T_Q::DOF;

	std::vector<Eigen::Matrix<T_SCALAR,DOF,1>> lambda;
	std::vector<Eigen::Matrix<T_SCALAR,DOF,1>> mu;

	bool m_initialized;

public:
	StepInternals (Abstract::LieProblem<T_SCALAR,T_Q,T_VEL>& problem, Epetra_Comm& comm)
	:	m_problem(&problem),
		m_comm(&comm)
	{
		int i,j;
		
		int myPID = m_comm->MyPID();
		int numProcs = m_comm->NumProc();

		size_t n_space_steps = m_problem->size(1);
		const int numGlobalEntries = DOF * n_space_steps;

		const int indexBase = 0;

		m_numLocalSystems = (n_space_steps+numProcs-myPID-1)/numProcs;
		int numLocalEntries = m_numLocalSystems*DOF;
		int* gblIndList = new int [numLocalEntries];

		for (i=0; i<m_numLocalSystems; i++) {
			for (j=0; j<DOF; j++) {
				gblIndList[i*DOF+j] = myPID*DOF + i*DOF*numProcs + j;
			}
		}

		m_map = new Epetra_Map(numGlobalEntries,numLocalEntries,gblIndList,indexBase,comm);

		if (gblIndList != NULL) {
			delete [] gblIndList;
			gblIndList = NULL;
		}

		m_initialGuess = Teuchos::rcp(new Epetra_Vector(*m_map));
		m_initialGuess->PutScalar(0.0);

		m_time_index = 0;

		lambda.resize(n_space_steps-1);
		mu.resize(n_space_steps-1);

		m_initialized = false;
	}

	~StepInternals ()
	{
		delete m_map;
	}

	Teuchos::RCP<Epetra_Vector>
	getInitialGuess ()
	{
		if (!m_initialized)
			return m_initialGuess;

		int i,j;
		int* gblIndList = m_map->MyGlobalElements();
		int numLocalEntries = m_map->NumMyElements();
		int space_index, time_index;
		if (numLocalEntries != DOF*m_numLocalSystems)
			std::cout << "Error: numLocalEntries != dof*m_numLocalSystems, " << numLocalEntries << " != " << DOF << "*" << m_numLocalSystems << std::endl;

		T_VEL vel;

		for (j=0; j<m_numLocalSystems; j++) {
			space_index = gblIndList[j*DOF]/DOF;
			if (m_time_index>0)
				time_index = m_time_index-1;
			else
				time_index = 0;
			vel = m_problem->vel(time_index,space_index,0);

			for (i=0; i<DOF; i++) {
				(*m_initialGuess)[j*DOF+i] = vel[i];
			}
		}

		return m_initialGuess;
	}

	bool
	initLimit ( )
	{
		// TODO: parallel
		// On suppose que la configuration du systeme est donné pour t = 0,h
		int j,k;
		double h;
		T_VEL xi;

		for(j=0; j<m_problem->size(1); j++) {
			xi = computeSpeed(0,j,0);
			computeSpeed(0,j,1);
			h = m_problem->step_size(0,j,0);
			//xi = (1.0/h)*T_VEL::cay_inv(m_problem->pos(0,j).inverse()*m_problem->pos(1,j));
			//m_problem->vel(0,j,0,xi);
			if (j<m_problem->size(1)-1)
				mu[j] = (h*xi).dCayRInv().transpose()*m_problem->dLdv0(xi);
			/*
			std::cout << "h: " << h << std::endl;
			std::cout << "fkk+1: " << std::endl << (m_problem->pos(0,j).inverse()*m_problem->pos(1,j)).matrix() << std::endl;
			std::cout << "xi:" << std::endl << xi << std::endl;
			std::cout << "mu [" << j << "]: " << std::endl << mu[j] << std::endl;
			*/
		}

		m_initialized = true;

		return true;
	}

	T_VEL
	computeSpeed (int i, int j, int dir)
	{
		/* ATTENTION: retourne la vitesse et l'ecrit*/
		// on suppose i in [0,n-1], j in [0,m-1], dir in [0,1]
		// on suppose connaitre la position en k et k+1 pour calculer la vitesse en k
		
		T_VEL v;
		double h;

		if (dir==0) {
			if (i==m_problem->size(0)-1) 
				v = T_VEL::Zero();
			else {
				h = m_problem->step_size(i,j,dir);
				v = (1.0/h)*T_VEL::cay_inv(m_problem->pos(i,j).inverse()*m_problem->pos(i,j+1));
			}
		}
		else {
			if (j==m_problem->size(1)-1)
				v = T_VEL::Zero();
			else {
				h = m_problem->step_size(i,j,dir);
				v = (1.0/h)*T_VEL::cay_inv(m_problem->pos(i,j).inverse()*m_problem->pos(i,j+1));
			}
		}

		m_problem->vel(i,j,dir,v);

		return v;
	}

	T_Q
	computePos (int i, int j, int dir)
	{
		// on suppose connaitre la vitesse et position en k pour calculer celle en k+1

		T_Q q;
		double h;
		T_VEL v;

		h = m_problem->step_size(i,j,dir);
		v = m_problem->vel(i,j,dir);
		q = m_problem->pos(i,j)*(h*v).cay();

		if (dir==0)
			m_problem->pos(i+1,j,q);
		else
			m_problem->pos(i,j+1,q);
	
		return q;
	}

	void
	setSolution (const Epetra_Vector& x)
	{
		int i,j;
		int* gblIndList = m_map->MyGlobalElements();
		int time_index = m_time_index, space_index;
		T_VEL xi;
		Eigen::Matrix<double,DOF,1> tmp_xi;

		for (j=0; j<m_numLocalSystems; j++) {
			space_index = gblIndList[j*DOF]/DOF;
			for (i=0; i<DOF; i++)
				tmp_xi[i] = x[i+j*DOF];
			xi = T_VEL(tmp_xi);
			m_problem->vel(time_index,space_index,0,xi);
			computePos(time_index,space_index,0);
		}
	}

	bool
	setupStep ( )
	{
		// TODO: parallel
		
		m_time_index++;
		int i = m_time_index, j;
		double h = m_problem->step_size(i,0,0);
		double s;

		Eigen::Matrix<double,6,1> mu_prev;
		T_VEL eta, xi_prev;

		for (j=0; j<m_problem->size(1)-1; j++) {
			s = m_problem->step_size(i,j,1);

			//eta =		(1.0/s)*T_VEL::cay_inv(m_problem->pos(i,j).inverse()*m_problem->pos(i,j+1));
			eta = computeSpeed(i,j,1);
			//xi_prev =	(1.0/h)*T_VEL::cay_inv(m_problem->pos(i-1,j).inverse()*m_problem->pos(i,j));
			lambda[j] = (s*eta).dCayRInv().transpose()*m_problem->dLdv1(eta);
			mu_prev =	mu[j];
			mu[j] =		T_VEL(mu_prev).Ad_star(m_problem->pos(i-1,j).inverse()*m_problem->pos(i,j)).toVector() + (h/s)*lambda[j];
			if (j>0)
				mu[j] += (-(h/s)*T_VEL(lambda[j-1]).Ad_star(m_problem->pos(i,j-1).inverse()*m_problem->pos(i,j))).toVector();
			//std::cout << "mu[" << j << "]:" << std::endl << mu[j] << std::endl;
		}

		return true;
	}

	bool
	computeF (const Epetra_Vector& x, Epetra_Vector& f, const FillType flag = Residual)
	{
		int i,j;
		int* gblIndList = m_map->MyGlobalElements();

		Eigen::Matrix<double,DOF,1> tmp_f;
		Eigen::Matrix<double,DOF,1> tmp_x;
		T_VEL xi;

		int time_index = m_time_index, space_index;
		double s, h = m_problem->step_size(time_index,0,0);

		for (j=0; j<m_numLocalSystems; j++) {
			space_index = gblIndList[j*DOF]/DOF;

			for (i=0; i<DOF; i++)
				tmp_x[i] = x[i+j*DOF];
			xi = T_VEL(tmp_x);

			//std::cout << "Asking for index " << space_index << " in mu of size " << mu.size() << std::endl;
			if (space_index<m_problem->size(1)-1)
				tmp_f =	(h*xi).dCayRInv().transpose()*m_problem->dLdv0(xi) - mu[space_index];
			else
				tmp_f = Eigen::Matrix<double,DOF,1>::Zero();

			//tmp_f = Eigen::Matrix<double,DOF,1>::Zero();

			for (i=0; i<DOF; i++) {
				f[j*DOF+i] = tmp_f[i];
			}
		}

		return true;
	}

};

} // namespace BiVariational

#endif
