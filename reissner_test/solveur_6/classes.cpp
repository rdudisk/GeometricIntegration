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


/* RigidBody */

Eigen::Matrix<double,6,6>& RigidBody::Inertia () { return m_Inertia; }
void RigidBody::Inertia (Eigen::Matrix<double,6,6> val) { m_Inertia = val; }
Eigen::Matrix<double,6,6>& RigidBody::Constraint () { return m_Constraint; }
void RigidBody::Constraint (Eigen::Matrix<double,6,6> val) { m_Constraint = val; }

void RigidBody::setInertia (double area, double rho)
{	// la formule n'est peut-être pas exacte, mais donne un ordre de grandeur
	double lm = area*rho; // masse lineaire
	Vec6 v;
	v << 0.5*lm*area, 0.25*lm*area, 0.25*lm*area, lm, lm, lm;
	m_Inertia = v.asDiagonal();
}

void RigidBody::setConstraint (double area, double young, double poisson)
{
	double G = young/(2.0+2.0*poisson);
	Vec6 v;
	v << G*(m_Inertia[1]+m_Inertia[2]), young*m_Inertia[1], young*m_Inertia[2], 
		young*area, G*area, G*area;
	m_Constraint = v.asDiagonal();
}

double RigidBody::coeffCFL (double young, double poisson, double rho, double alpha)
{ return 1.0/(alpha*sqrt(young*(1.0-poisson)/(rho*(1.0+poisson)*(1.0-2.0*poisson)))); }

void
RigidBody::writeCSVFile (const std::string filename, bool header)
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
		of << "#i,j,t,s,x,y,z,u1,u2,u3,v1,v2,v3,k,b,m1,m2,m3,m4,m5,m6" << std::endl;

	for (i=0; i<n_time-1; i++) {
		T = i*h;
		for (j=0; j<n_space; j++) {
			S = j*l;
			p = this->pos(i,j);
			x = p.trans();
			u = p.rotateVector(E2);
			v = p.rotateVector(E3);
			v_xi = this->vel_time(i,j).toVector();
			// osef si j=n_space-1
			v_eps = this->vel_space(i,j).toVector();
			kinetic = 0.5*(v_xi.dot((this->Inertia())*v_xi));
			if (j==n_space-1) bending = 0.0;
			else bending = 0.5*((v_eps-E4).dot((this->Constraint())*(v_eps-E4)));
			momentum = l*Algebra::static_Ad_star(p.inverse(),Algebra(this->mom_time(i,j)));
			/*
			c = xi.trans();
			d = eta.trans();
			if (j<n_space-1) eta = this->vel(i,j,1);
			else eta = Algebra::Zero();
			if (i<n_time-1) xi = this->vel(i,j,0);
			else xi = Algebra::Zero();
			*/
			of	<< i << "," << j << ","
				<< T << "," << S << ","
				<< x[0] << "," << x[1] << "," << x[2] << ","
				<< u[0] << "," << u[1] << "," << u[2] << ","
				<< v[0] << "," << v[1] << "," << v[2];
				//<< "," << kinetic << "," << bending;
			for (int k=0; k<6; k++)
				of << "," << momentum[k];
			of << std::endl;
				//<< c[0] << "," << c[1] << "," << c[2] << ","
				//<< d[0] << "," << d[1] << "," << d[2] << std::endl;
		}
	}

	of.close();
}

/* SolveMe */

void
SolveMe::setData (double h, Vec6 mu)
{ m_h = h; m_mu = mu; }

void
SolveMe::setx0 (Vec6 x0)
{ m_x0 = x0; }

double SolveMe::h () const { return m_h; }
Vec6 SolveMe::mu () const { return m_mu; }

const NOXVector<6>
SolveMe::getInitialGuess ()
{
	NOXVector<6> ret(m_x0);
	return ret;
}

bool
SolveMe::computeF (NOXVector<6>& f, const NOXVector<6>& x)
{
	Eigen::Matrix<double,6, 6> J;
	J = this->m_problem.Inertia();

	// TODO: adapter
	double lambda = 1.0+pow(OM.norm(),2);
	Eigen::Matrix<double,3,3> OM_hat = SO3::Algebra<double>(OM).toRotationMatrix();
	Eigen::Matrix<double,3,3> A = Eigen::Matrix<double,3,3>::Identity()+OM_hat;
	Eigen::Matrix<double,3,3> Ainv =
		(1.0/lambda)*(Eigen::Matrix<double,3,3>::Identity()-OM_hat+OM*OM.transpose());
	Vec3 Gamma = J2.inverse()*Ainv*M2;
	f = lambda*Ainv.transpose()*J1*OM + A*SO3::Algebra<double>(Gamma).toRotationMatrix()*Ainv*M2 - M1;

	return true;
}

bool
SolveMe::computeJacobian (Eigen::Matrix<double,6,6>& J, const NOXVector<6>& x)
{
	// TODO: adapter
	// voir cahier 5 p. 59
	double nOM2 = pow(OM.norm(),2);
	double lambda = 1.0+nOM2;
	Eigen::Matrix<double,3,3> J1, J2, A, B, Binv, C, Q, I, OM_hat, dGdOM;
	Eigen::Matrix<double,3,1> V, Gamma;
	J1 = this->m_problem.Inertia().block(0,0,3,3);
	J2 = this->m_problem.Inertia().block(3,0,3,3);
	OM_hat = SO3::Algebra<double>(OM).toRotationMatrix();
	I = Eigen::Matrix<double,3,3>::Identity();
	A = I+OM_hat+OM*OM.transpose();
	B = I+OM_hat;
	Binv = (1.0/lambda)*(I-OM_hat+OM*OM.transpose());
	C = SO3::Algebra<double>(M2).toRotationMatrix() + OM*M2.transpose() + OM.transpose()*M2*I;
	Gamma = J2.inverse()*Binv*M2;
	V = A.transpose()*M2;
	dGdOM = (1.0/lambda)*J2.inverse()*(C-(1.0/(nOM2*lambda))*V*(OM.transpose()));
	for (int j=0; j<3; j++)
		Q.col(j) = SO3::Algebra<double>(dGdOM.col(j)).toRotationMatrix()*V;
	J = A*J1 - SO3::Algebra<double>(J1*OM).toRotationMatrix()
		+ OM*((J1*OM).transpose()) + ((OM.transpose())*J1*OM)*I
		+ (1.0/lambda)*( -SO3::Algebra<double>(SO3::Algebra<double>(Gamma).toRotationMatrix()*(A.transpose())*M2).toRotationMatrix()
				+B*(Q + SO3::Algebra<double>(Gamma).toRotationMatrix()*((-1.0/(nOM2*lambda))*(A.transpose()*M2*(OM.transpose())) + C)));
	return true;
}
