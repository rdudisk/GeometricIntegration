#include <iostream>
#include <fstream>
#include <cmath>

#include "header.hpp"

/* InverseLegendre */

void
InverseLegendre::setData (double h, Vec6 mu)
{ m_h = h; m_M = 0.5*h*mu; }

const NOXVector<6>
InverseLegendre::getInitialGuess ()
{
	NOXVector<6> ret = this->m_problem.Inertia().inverse()*this->m_M;
	return ret;
}

bool
InverseLegendre::computeF (NOXVector<6>& f, const NOXVector<6>& X)
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
InverseLegendre::computeJacobian (Eigen::Matrix<double,6,6>& J, const NOXVector<6>& X)
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

bool
InverseLegendre::computeSolution () {
	bool success, verbose;
	verbose = false;
	try {
		m_solver->reset(this->getInitialGuess());
		NOX::StatusTest::StatusType status = this->m_solver->solve();
		const NOXGroup<NOXVector<6>,1>& solnGrp = dynamic_cast<const NOXGroup<NOXVector<6>,1>&>(this->m_solver->getSolutionGroup());
		const NOXVector<6>& c_x = dynamic_cast<const NOXVector<6>&>(solnGrp.getX());

		if (status == NOX::StatusTest::Failed || status == NOX::StatusTest::Unconverged)
			success = false;
		else success = true;

		this->m_xsol = Algebra((2.0/m_h)*c_x);
	} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
	if (!success) {
		std::cout << "Fail with" << std::endl;
		std::cout << "M:  " << m_M << std::endl;
		std::cout << "X:  " << 0.5*m_h*m_xsol.vector() << std::endl;
	}
	return success;
}

Algebra
InverseLegendre::getAlgebraSolution () {
	return this->m_xsol;
}

/* RestrictedInverseLegendre */

void
RestrictedInverseLegendre::setData (double h, Vec6 mu)
{ m_h = h; m_mu = mu; }

const NOXVector<2>
RestrictedInverseLegendre::getInitialGuess ()
{
	NOXVector<2> ret;
	ret[0] = m_mu[1]/m_problem.Inertia()(1,1);
	ret[1] = m_mu[2]/m_problem.Inertia()(2,2);
	return ret;
}

bool
RestrictedInverseLegendre::computeF (NOXVector<2>& f, const NOXVector<2>& x)
{
	// voir cahier 6 p. 170
	double j2,j3,w2,w3,m2,m3,h2,h_coeff;
	j2     = this->m_problem.Inertia()(1,1);
	j3     = this->m_problem.Inertia()(2,2);
	w2     = x[0];
	w3     = x[1];
	m2     = m_mu[1];
	m3     = m_mu[2];
	h_coeff = 0.25*m_h*m_h;

	f[0] = j2*w2*(1.0+h_coeff*w2*(w2+w3)) - m2;
	f[1] = j3*w3*(1.0+h_coeff*w3*(w2+w3)) - m3;
	
	return true;
}

bool
RestrictedInverseLegendre::computeJacobian (Eigen::Matrix<double,2,2>& J, const NOXVector<2>& x)
{
	// voir cahier 6 p. 171
	double j2,j3,w2,w3,m2,m3,h2,h_coeff;
	j2     = this->m_problem.Inertia()(1,1);
	j3     = this->m_problem.Inertia()(2,2);
	w2     = x[0];
	w3     = x[1];
	m2     = m_mu[1];
	m3     = m_mu[2];
	h_coeff = 0.25*m_h*m_h;

	J(0,0) = j2*(1.0+h_coeff*w2*(2.0*w3+3.0*w2));
	J(0,1) = j2*h_coeff*w2*w2;
	J(1,0) = j3*h_coeff*w3*w3;
	J(1,1) = j3*(1.0+h_coeff*w3*(2.0*w2+3.0*w3));

	return true;
}

bool
RestrictedInverseLegendre::computeSolution () {
	bool success, verbose;
	verbose = false;
	try {
		m_solver->reset(this->getInitialGuess());
		NOX::StatusTest::StatusType status = this->m_solver->solve();
		const NOXGroup<NOXVector<2>,1>& solnGrp = dynamic_cast<const NOXGroup<NOXVector<2>,1>&>(this->m_solver->getSolutionGroup());
		const NOXVector<2>& c_x = dynamic_cast<const NOXVector<2>&>(solnGrp.getX());

		if (status == NOX::StatusTest::Failed || status == NOX::StatusTest::Unconverged)
			success = false;
		else success = true;

		Eigen::Matrix<double,6,1> v;
		v << 0.0, c_x[0], c_x[1], 0.0, 0.0, 0.0;

		this->m_xsol = Algebra(v);
	} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
	if (!success) {
		std::cout << "Failed restricted inverse legendre!" << std::endl;
	}
	return success;
}

Algebra
RestrictedInverseLegendre::getAlgebraSolution () {
	return this->m_xsol;
}

