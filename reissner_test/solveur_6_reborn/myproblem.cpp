#include "header.hpp"

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

void RigidBody::setTensors (double radius, double rho, double young, double shear)
{
	double area = RigidBody::compute_area(radius);
	double lm   = area*rho; // masse lineique
	double Irho = 0.5*area*radius*radius;
	double Ia   = 0.5*Irho;
	Vec6 v;

	v << rho*Irho, rho*Ia, rho*Ia, lm, lm, lm;
	m_Inertia = v.asDiagonal();

	v<< shear*Irho, young*Ia, young*Ia, young*area, shear*area, shear*area;
	m_Constraint = v.asDiagonal();
}

double RigidBody::coeffCFL (double young, double shear, double rho, double alpha)
{
	double cfl, lambda;

	lambda = shear*(young-2.0*shear)/(3.0*shear-young);
	cfl	   = alpha/sqrt((lambda+2.0*shear)/rho);

	return cfl;
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

