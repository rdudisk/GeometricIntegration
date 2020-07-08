#ifndef DEF_BIVARIATIONAL_NORMAL_SOLVEME
#define DEF_BIVARIATIONAL_NORMAL_SOLVEME

namespace BiVariational {

template <typename T_SCALAR,
		 typename T_Q,
		 typename T_VEL>
class NormalStepInternals : public ::Abstract::NOXStep<T_VEL,1>
{
	//using MyEpetraOperator::m_jacobian;

private:
	Abstract::Problem<T_SCALAR,T_Q,T_VEL>* m_problem;
	Epetra_Map* m_map;
	Epetra_Comm* m_comm;

	Teuchos::RCP<Epetra_Vector> m_initialGuess;
	//Teuchos::RCP<Epetra_VbrMatrix> m_jacobianMatrix;
	Teuchos::RCP<Epetra_VbrMatrix> m_jacobian;

	int m_numLocalSystems;
	int m_time_index;

	static const int DOF = T_Q::DOF;

	std::vector<Eigen::Matrix<T_SCALAR,DOF,1>> lambda;
	std::vector<Eigen::Matrix<T_SCALAR,DOF,1>> mu;

	bool m_initialized;

public:
	NormalStepInternals (Abstract::Problem<T_SCALAR,T_Q,T_VEL>& problem, Epetra_Comm& comm)
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

		m_initialGuess = Teuchos::rcp(new Epetra_Vector(*m_map));
		m_initialGuess->PutScalar(0.0);

		// TODO Matrice: vérifier
		int* gblBlockIndList = new int [m_numLocalSystems];
		for (i=0; i<m_numLocalSystems; i++)
			gblBlockIndList[i] = myPID + i*numProcs;
		int* rowSizeList = new int [m_numLocalSystems];
		for (i=0; i<m_numLocalSystems; i++)
			rowSizeList[i] = DOF;
		Epetra_BlockMap blockMap(n_space_steps,m_numLocalSystems,gblBlockIndList,rowSizeList,indexBase,comm);

		m_jacobian = Teuchos::rcp(new Epetra_VbrMatrix(Copy,blockMap,1));
		//m_jacobian = new Epetra_VbrMatrix(Copy,blockMap,1);

		int globalNode;
		int* indices = new int [1];
		for (i=0; i<m_numLocalSystems; i++) {
			globalNode = gblBlockIndList[i];
			indices[0] = globalNode;
			m_jacobian->BeginInsertGlobalValues(globalNode,1,indices);
			double* values = new double [DOF*DOF];
			for (j=0; j<DOF*DOF; j++) {
				values[j] = 1.0*j;
			}
			m_jacobian->SubmitBlockEntry(values,6,6,6);
			m_jacobian->EndSubmitEntries();
		}

		m_jacobian->FillComplete();

		m_time_index = 0;

		lambda.resize(n_space_steps-1);
		mu.resize(n_space_steps-1);

		m_initialized = false;

		delete [] gblIndList;
		delete [] gblBlockIndList;
		delete [] rowSizeList;
		delete [] indices;
	}

	~NormalStepInternals ()
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

	bool computeJacobian (const Epetra_Vector& x, Epetra_Operator& Jac)
	{
		Epetra_VbrMatrix* tmp_matrix = dynamic_cast<Epetra_VbrMatrix*>(&Jac);
		//m_jacobian =
		int i,j,k,l;
		int* gblIndList = m_map->MyGlobalElements();
		int myPID = m_comm->MyPID();
		int numProcs = m_comm->NumProc();

		Eigen::Matrix<double,DOF,DOF> tmp_jac, A, jac_problem;
		Eigen::Matrix<double,DOF,1> tmp_x;
		T_VEL xi;

		Eigen::Matrix<double,3,1>
			E1 = SO3::Algebra<double>::GeneratorVector(0),
			E2 = SO3::Algebra<double>::GeneratorVector(1),
			E3 = SO3::Algebra<double>::GeneratorVector(2);
		std::vector<Eigen::Matrix<double,3,1>> v_E;
		v_E.push_back(E1); v_E.push_back(E2); v_E.push_back(E3);
		Eigen::Matrix<double,3,3>
			hat_E1 = SO3::Algebra<double>::GeneratorMatrix(0),
			hat_E2 = SO3::Algebra<double>::GeneratorMatrix(1),
			hat_E3 = SO3::Algebra<double>::GeneratorMatrix(2);
		std::vector<Eigen::Matrix<double,3,3>> v_hat_E;
		v_hat_E.push_back(hat_E1); v_hat_E.push_back(hat_E2); v_hat_E.push_back(hat_E3);

		int time_index = m_time_index, space_index;
		double s, h = m_problem->step_size(time_index,0,0);

		int globalNode;
		int* gblBlockIndList = new int [m_numLocalSystems];
		for (i=0; i<m_numLocalSystems; i++)
			gblBlockIndList[i] = myPID + i*numProcs;
		int* indices = new int [1];
		double* values = new double [DOF*DOF];

		for (j=0; j<m_numLocalSystems; j++) {
			space_index = gblIndList[j*DOF]/DOF;

			for (i=0; i<DOF; i++)
				tmp_x[i] = x[i+j*DOF];
			xi = T_VEL(tmp_x);

			if (space_index<m_problem->size(1)-1) {
				jac_problem = m_problem->d2Ldv0(xi);
				for (i=0; i<3; i++) {
					// TODO: vérifier indices matrice
					A.block(0,0,3,3) = (-h/2.0)*v_hat_E[i]+(h*h/2.0)*(v_E[i]*xi.rot().transpose()+xi.rot()*v_E[i].transpose());
					A.block(0,3,3,3) = Eigen::Matrix<double,3,3>::Zero();
					A.block(3,0,3,3) = (h*h/4.0)*v_hat_E[i]*SO3::Algebra<double>(xi.trans()).toRotationMatrix();
					A.block(3,3,3,3) = (-h/2.0)*v_hat_E[i];
					tmp_jac.col(i) = A.transpose()*m_problem->dLdv0(xi) + (h*xi).dCayRInv().transpose()*jac_problem.col(i);
				}
				for (i=3; i<DOF; i++) {
					// TODO: vérifier indices matrice
					A.block(0,0,3,6) = Eigen::Matrix<double,3,6>::Zero();
					A.block(3,3,3,3) = Eigen::Matrix<double,3,3>::Zero();
					A.block(3,0,3,3) = (-h/2.0)*(Eigen::Matrix<double,3,3>::Identity()-(h/2.0)*xi.rotationMatrix())*v_hat_E[i-3];
					tmp_jac.col(i) = A.transpose()*m_problem->dLdv0(xi) + (h*xi).dCayRInv().transpose()*jac_problem.col(i);
				}
			}
			else {
					tmp_jac = Eigen::Matrix<double,6,6>::Identity();
			}

			/* Fill the actual matrix */
			for (i=0; i<m_numLocalSystems; i++) {
				globalNode = gblBlockIndList[i];
				indices[0] = globalNode;
				tmp_matrix->BeginReplaceGlobalValues(globalNode,1,indices);
				values = tmp_jac.data();
				tmp_matrix->SubmitBlockEntry(values,6,6,6);
				tmp_matrix->EndSubmitEntries();
			}

			tmp_matrix->FillComplete();
		}

		return true;
	}

	Teuchos::RCP<Epetra_VbrMatrix>
	jacobian () const
	{ return m_jacobian; }

};

} // namespace BiVariational

#endif
