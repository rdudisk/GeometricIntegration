#include <iostream>
#include <fstream>
#include <cmath>

#include "classes.hpp"

/* MultiSyst */

M MultiSyst::base_time (void) const { return m_base_time; }
void MultiSyst::base_time(M b) { m_base_time = b; }
M MultiSyst::base_space (void) const { return m_base_space; }
void MultiSyst::base_space(M b) { m_base_space = b; }
Group MultiSyst::pos ( ) const { return m_pos; }
void MultiSyst::pos (Group p) { m_pos = p; }
Algebra MultiSyst::vel_time ( ) const { return m_vel_time; }
void MultiSyst::vel_time (Algebra v) { m_vel_time = v; }
Algebra MultiSyst::vel_space ( ) const { return m_vel_space; }
void MultiSyst::vel_space (Algebra v) { m_vel_space = v; }
Vec6 MultiSyst::mom_time ( ) const { return m_mom_time; }
void MultiSyst::mom_time (Vec6 m) { m_mom_time = m; }
Vec6 MultiSyst::mom_space ( ) const { return m_mom_space; }
void MultiSyst::mom_space (Vec6 m) { m_mom_space = m; }

/* DiscMultiSyst */

void
DiscMultiSyst::setSize (const size_t i, const size_t j)
{
	m_size[0] = i;
	m_size[1] = j;
	m_node.resize(i*j);

	MultiSyst default_node(0.0,0.0);
	for (int k=0; k<i*j; k++) m_node[k] = default_node;
}

MultiSyst const& DiscMultiSyst::operator[] (size_t index) const {
	if (index<0 || index>=m_node.size())
		throw std::out_of_range("Error : DiscMultiSyst[] index out of range");
	return m_node[index];
}

MultiSyst& DiscMultiSyst::operator[] (size_t index) {
	if (index<0 || index>=m_node.size())
		throw std::out_of_range("Error : DiscMultiSyst[] index out of range");
	return m_node[index];
}

size_t DiscMultiSyst::size(const size_t dir) const { return m_size[dir]; }

const size_t DiscMultiSyst::getIndex (const size_t i, const size_t j) const
{ return i*m_size[1]+j; }

M DiscMultiSyst::step_size (const size_t dir) const
{
	M h = base_time(1)-base_time(0);
	if (dir==1) { h = base_space(1)-base_space(0); }
	return h;
}

M DiscMultiSyst::base_time (const size_t& i) const {return m_node[getIndex(i,0)].base_time(); }
M DiscMultiSyst::base_space (const size_t& j) const {return m_node[getIndex(0,j)].base_space(); }

Group DiscMultiSyst::pos (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].pos(); }
void DiscMultiSyst::pos (const size_t& i, const size_t& j, const Group& p) { m_node[getIndex(i,j)].pos(p); }

Algebra DiscMultiSyst::vel_time (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].vel_time(); }
void DiscMultiSyst::vel_time (const size_t& i, const size_t& j, const Algebra& v) { m_node[getIndex(i,j)].vel_time(v); }

Algebra DiscMultiSyst::vel_space (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].vel_space(); }
void DiscMultiSyst::vel_space (const size_t& i, const size_t& j, const Algebra& v) { m_node[getIndex(i,j)].vel_space(v); }

Vec6 DiscMultiSyst::mom_time (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].mom_time(); }
void DiscMultiSyst::mom_time (const size_t& i, const size_t& j, const Vec6& m) { m_node[getIndex(i,j)].mom_time(m); }

Vec6 DiscMultiSyst::mom_space (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].mom_space(); }
void DiscMultiSyst::mom_space (const size_t& i, const size_t& j, const Vec6& m) { m_node[getIndex(i,j)].mom_space(m); }

const unsigned int DiscMultiSyst::dof ( ) { return Group::DOF; }

void DiscMultiSyst::baselinstep (M t_inf_lim, M t_step_size, M s_inf_lim, M s_step_size)
{
	m_node = std::vector<MultiSyst>();
	for (int i=0; i<m_size[0]; i++) { for (int j=0; j<m_size[1]; j++) {
			m_node.push_back(MultiSyst(t_inf_lim+i*t_step_size,s_inf_lim+j*s_step_size));
	} }
}


/* SparseDiscMultiSyst */

void
SparseDiscMultiSyst::setSize (const size_t j)
{
	// m_size[0] = TIME_SIZE;
	m_size[1] = j;
	m_node.resize(m_size[0]*m_size[1]);

	MultiSyst default_node(0.0,0.0);
	for (int k=0; k<m_size[0]*m_size[1]; k++) m_node[k] = default_node;
}

void
SparseDiscMultiSyst::setStepSize (const M h, const M l)
{
	m_step_size[0] = h;
	m_step_size[1] = l;
}

MultiSyst const& SparseDiscMultiSyst::operator[] (size_t index) const {
	if (index<0 || index>=m_size[0]*m_size[1])
		throw std::out_of_range("Error : DiscMultiSyst[] index out of range");
	return m_node[index];
}

MultiSyst& SparseDiscMultiSyst::operator[] (size_t index) {
	if (index<0 || index>=m_size[0]*m_size[1])
		throw std::out_of_range("Error : DiscMultiSyst[] index out of range");
	return m_node[index];
}

size_t SparseDiscMultiSyst::size(const size_t dir) const { return m_size[dir]; }

const size_t SparseDiscMultiSyst::getIndex (const size_t i, const size_t j) const
{ return (i%TIME_SIZE)*m_size[1]+j; }

/*
M SparseDiscMultiSyst::step_size (const size_t dir) const
{
	M h = base_time(1)-base_time(0);
	if (dir==1) { h = base_space(1)-base_space(0); }
	return h;
}*/

M SparseDiscMultiSyst::base_time (const size_t& i) const {return m_node[getIndex(i,0)].base_time(); }
M SparseDiscMultiSyst::base_space (const size_t& j) const {return m_node[getIndex(0,j)].base_space(); }

Group SparseDiscMultiSyst::pos (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].pos(); }
void SparseDiscMultiSyst::pos (const size_t& i, const size_t& j, const Group& p) { m_node[getIndex(i,j)].pos(p); }

Algebra SparseDiscMultiSyst::vel_time (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].vel_time(); }
void SparseDiscMultiSyst::vel_time (const size_t& i, const size_t& j, const Algebra& v) { m_node[getIndex(i,j)].vel_time(v); }

Algebra SparseDiscMultiSyst::vel_space (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].vel_space(); }
void SparseDiscMultiSyst::vel_space (const size_t& i, const size_t& j, const Algebra& v) { m_node[getIndex(i,j)].vel_space(v); }

Vec6 SparseDiscMultiSyst::mom_time (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].mom_time(); }
void SparseDiscMultiSyst::mom_time (const size_t& i, const size_t& j, const Vec6& m) { m_node[getIndex(i,j)].mom_time(m); }

Vec6 SparseDiscMultiSyst::mom_space (const size_t& i, const size_t& j) const { return m_node[getIndex(i,j)].mom_space(); }
void SparseDiscMultiSyst::mom_space (const size_t& i, const size_t& j, const Vec6& m) { m_node[getIndex(i,j)].mom_space(m); }

const unsigned int SparseDiscMultiSyst::dof ( ) { return Group::DOF; }

void SparseDiscMultiSyst::setCSV (const std::string filename)
{
	this->csv_file.open(filename,std::ios_base::trunc);
	this->csv_file << "i,j,t,s,x,y,z,u1,u2,u3,v1,v2,v3,k,b,m1,m2,m3,m4,m5,m6" << std::endl;
}

void SparseDiscMultiSyst::endCSV ()
{
	this->csv_file.close();
}


/* RigidBody */

double RigidBody::compute_area (double radius)
{ return M_PI*radius*radius; }

Eigen::Matrix<double,6,6>& RigidBody::Inertia () { return m_Inertia; }
void RigidBody::Inertia (Eigen::Matrix<double,6,6> val) { m_Inertia = val; }
Eigen::Matrix<double,6,6>& RigidBody::Constraint () { return m_Constraint; }
void RigidBody::Constraint (Eigen::Matrix<double,6,6> val) { m_Constraint = val; }

void RigidBody::setInertia (double radius, double rho)
{
	double area = RigidBody::compute_area(radius);
	double lm = area*rho; // masse lineique
	double Ia = 0.5*lm*area;
	Vec6 v;
	v << Ia, 0.5*Ia, 0.5*Ia, lm, lm, lm;
	m_Inertia = v.asDiagonal();
}

void RigidBody::setConstraint (double radius, double young, double poisson)
{
	double area = RigidBody::compute_area(radius);
	double G = young/(2.0+2.0*poisson);
	Vec6 v;
	v<< G*m_Inertia(0,0), young*m_Inertia(1,1), young*m_Inertia(2,2), 
		young*area, G*area, G*area;
	m_Constraint = v.asDiagonal();
}

double RigidBody::coeffCFL (double young, double poisson, double rho, double alpha)
{
	double cfl, cfl_factor, lame_lambda, lame_mu;

	lame_lambda = young*poisson/((1.0+poisson)*(1.0-2.0*poisson));
	lame_mu     = young/(2.0*(1.0+poisson));
	cfl	        = alpha/sqrt((lame_lambda+2.0*lame_mu)/rho);

	return cfl;
	//return 1.0/(alpha*sqrt(young*(1.0-poisson)/(rho*(1.0+poisson)*(1.0-2.0*poisson)))); }
}

void
RigidBody::updateCSV (int i, double w_resample)
{
	// Ecrit les données pour current_i%4 dans le fichier CSV
	// current_i      : l'indice i actuel
	// w_resample     : la periode entre deux echantillons
	// this->m_last_i : le dernier indice i pour lequel on a ecrit dans le fichier
	
	double h = this->m_step_size[0];
	double l = this->m_step_size[1];

	int Last_i, Current_i;
	// Voir dessin cahier 6. p. 118
	if (w_resample!=0.0 && m_last_i!=-1) {
		Last_i    = (int) (h*((double)m_last_i)/w_resample);
		Current_i = (int) (h*((double)i)/w_resample);
		if (Last_i==Current_i) {
			return;
		}
	}

	int n_space = this->m_size[1];
	double T,S;
	Algebra xi, eta;
	Group p;
	Eigen::Matrix<double,3,1> x,u,v,c,d,E2,E3;
	E2 << 0, l, 0;
	E3 << 0, 0, l;
	Eigen::Matrix<double,6,1> E4, v_xi, v_eps;
	E4 << 0,0,0,1,0,0;
	double kinetic, bending;
	Algebra momentum;

	T = i*h;
	for (int j=0; j<n_space; j++) {
		S = j*l;
		p = this->pos(i,j);
		x = p.translationVector();
		u = p.rotateVector(E2);
		v = p.rotateVector(E3);
		v_xi = this->vel_time(i,j).vector();
		// osef si j=n_space-1
		v_eps = this->vel_space(i,j).vector();
		kinetic = 0.5*(v_xi.dot((this->Inertia())*v_xi));
		if (j==n_space-1) bending = 0.0;
		else bending = 0.5*((v_eps-E4).dot((this->Constraint())*(v_eps-E4)));
		momentum = l*Algebra::static_Ad_star(p.inverse(),Algebra(this->mom_time(i,j)));
		this->csv_file
			<< i << "," << j << ","
			<< T << "," << S << ","
			<< x[0] << "," << x[1] << "," << x[2] << ","
			<< u[0] << "," << u[1] << "," << u[2] << ","
			<< v[0] << "," << v[1] << "," << v[2]
			<< "," << kinetic << "," << bending;
		for (int k=0; k<6; k++)
			this->csv_file << "," << momentum[k];
		this->csv_file << std::endl;
	}

	m_last_i = i;
}

/*
void
RigidBody::writeCSVFile (const std::string filename, bool header, int resample)
{
	std::ofstream of;
	of.open(filename,std::ios_base::trunc);

	int i,j;
	int n_time = this->size(0);
	int n_space = this->size(1);
	double T,S;
	double h = this->step_size(0);
	double l = this->step_size(1);
	Algebra xi, eta;
	Group p;
	Eigen::Matrix<double,3,1> x,u,v,c,d,E2,E3;
	E2 << 0, l, 0;
	E3 << 0, 0, l;
	Eigen::Matrix<double,6,1> E4, v_xi, v_eps;
	E4 << 0,0,0,1,0,0;
	double kinetic, bending;
	Algebra momentum;
	
	if (header)
		of << "i,j,t,s,x,y,z,u1,u2,u3,v1,v2,v3,k,b,m1,m2,m3,m4,m5,m6" << std::endl;

	int i_step;
	if (resample==0) 
		i_step = 1;
	else
		i_step = n_time/resample;

	for (i=0; i<n_time-1; i+=i_step) {
		T = i*h;
		for (j=0; j<n_space; j++) {
			S = j*l;
			p = this->pos(i,j);
			x = p.translationVector();
			u = p.rotateVector(E2);
			v = p.rotateVector(E3);
			v_xi = this->vel_time(i,j).vector();
			// osef si j=n_space-1
			v_eps = this->vel_space(i,j).vector();
			kinetic = 0.5*(v_xi.dot((this->Inertia())*v_xi));
			if (j==n_space-1) bending = 0.0;
			else bending = 0.5*((v_eps-E4).dot((this->Constraint())*(v_eps-E4)));
			momentum = l*Algebra::static_Ad_star(p.inverse(),Algebra(this->mom_time(i,j)));
			of	<< i << "," << j << ","
				<< T << "," << S << ","
				<< x[0] << "," << x[1] << "," << x[2] << ","
				<< u[0] << "," << u[1] << "," << u[2] << ","
				<< v[0] << "," << v[1] << "," << v[2]
				<< "," << kinetic << "," << bending;
			for (int k=0; k<6; k++)
				of << "," << momentum[k];
			of << std::endl;
				//<< c[0] << "," << c[1] << "," << c[2] << ","
				//<< d[0] << "," << d[1] << "," << d[2] << std::endl;
		}
	}

	of.close();
}
*/

/* SolveMe */

void
SolveMe::setData (double h, Vec6 mu, Vec6 xi)
{ m_h = h; m_M = 0.5*h*mu; m_X0 = 0.5*h*xi; }

//double SolveMe::h () const { return m_h; }
//Vec6 SolveMe::M () const { return m_M; }

const NOXVector<6>
SolveMe::getInitialGuess ()
{
	//NOXVector<3> ret((1.0+1.0/this->m_h)*this->m_OM1-(1.0/this->m_h)*this->m_OM0);
	//NOXVector<3> ret(2*m_OM1-m_OM0);
	//NOXVector<6> ret(m_X0);
	NOXVector<6> ret = this->m_problem.Inertia().inverse()*this->m_M;
	return ret;
}

bool
SolveMe::computeF (NOXVector<6>& f, const NOXVector<6>& X)
{
	// Achtung! X = (h/2)*xi
	// de meme pour M
	Eigen::Matrix<double,3,1> F1, F2, OM, G;
	Eigen::Matrix<double,3,3> Id, J1, J2, OM_hat, G_hat;
	OM     = X.head(3);
	G      = X.tail(3);
	Id     = Eigen::Matrix<double,3,3>::Identity();
	J1     = this->m_problem.Inertia().block(0,0,3,3);
	J2     = this->m_problem.Inertia().block(3,3,3,3);
	OM_hat = SO3::Algebra<double>(OM).rotationMatrix();
	G_hat  = SO3::Algebra<double>(G).rotationMatrix();

	F1     = (Id+OM_hat+OM*OM.transpose())*J1*OM + (Id+OM_hat)*G_hat*J2*G;
	F2     = (Id+OM_hat)*J2*G;
	
	f.head(3) = (2.0/m_h)*F1;
	f.tail(3) = (2.0/m_h)*F2;
	f -= (2.0/m_h)*this->m_M;

	return true;
}

bool
SolveMe::computeJacobian (Eigen::Matrix<double,6,6>& J, const NOXVector<6>& X)
{
	// voir cahier 6 p. 110
	Eigen::Matrix<double,3,1> OM, G;
	Eigen::Matrix<double,3,3> JF1_OM, JF1_G, JF2_OM, JF2_G, Id, J1, J2, OM_hat, G_hat;
	OM     = X.head(3);
	G      = X.tail(3);
	Id     = Eigen::Matrix<double,3,3>::Identity();
	J1     = this->m_problem.Inertia().block(0,0,3,3);
	J2     = this->m_problem.Inertia().block(3,3,3,3);
	OM_hat = SO3::Algebra<double>(OM).rotationMatrix();
	G_hat  = SO3::Algebra<double>(G).rotationMatrix();
	
	JF1_OM  = (Id+OM_hat+OM*OM.transpose())*J1;
	JF1_OM -= SO3::Algebra<double>(J1*OM).rotationMatrix();
	JF1_OM += OM*(J1*OM).transpose();
	JF1_OM += (OM.transpose()*J1*OM)*Id;
	JF1_OM -= SO3::Algebra<double>(G_hat*J2*G).rotationMatrix();

	JF1_G   = (Id+OM_hat)*(G_hat*J2 -SO3::Algebra<double>(J2*G).rotationMatrix());

	JF2_OM  = (-1.0)*SO3::Algebra<double>(J2*G).rotationMatrix();

	JF2_G   = (Id+OM_hat)*J2;

	J.block(0,0,3,3) = (2.0/m_h)*JF1_OM;
	J.block(0,3,3,3) = (2.0/m_h)*JF1_G;
	J.block(3,0,3,3) = (2.0/m_h)*JF2_OM;
	J.block(3,3,3,3) = (2.0/m_h)*JF2_G;

	return true;
}

void
Displacement::flush ()
{
	for (int j=0; j<i; j++) {
		of << m_buffer[j] << std::endl;
	}
}

void
Displacement::insert (double y)
{
	m_buffer[i] = y;
	i++;
	if (i==BUFF_SIZE) {
		this->flush();
		i=0;
	}
}

void
Displacement::close ()
{
	this->of.close();
}
