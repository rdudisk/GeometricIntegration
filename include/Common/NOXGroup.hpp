#ifndef DEF_NOXGROUP
#define DEF_NOXGROUP

#include "NOX_Common.H"
#include "NOX_Abstract_Group.H"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "include/Common/NOXVector.hpp"
#include "include/Step/AbstractStep.hpp"

template <typename T_Q, int T_SIZE_MULTIVECT = 1>
class NOXGroup : public virtual NOX::Abstract::Group {
public:
	NOXGroup<T_Q,T_SIZE_MULTIVECT> (Abstract::Step<T_Q,T_SIZE_MULTIVECT>& step) :
		m_step(step),
		m_xVector(step.getInitialGuess())
	{
		resetIsValid();
	}

	~NOXGroup<T_Q,T_SIZE_MULTIVECT> ( )
	{ }

	NOX::Abstract::Group&
	operator= (const NOXGroup<T_Q,T_SIZE_MULTIVECT>& g)
	{
		if (this != &g) {
			m_xVector = g.m_xVector;

			m_isValidF = g.m_isValidF;
			m_isValidJacobian = g.m_isValidJacobian;
			m_isValidGradient = g.m_isValidGradient;
			m_isValidNewton = g.m_isValidNewton;

			if (m_isValidF)
				m_fVector = g.m_fVector;
			
			if (m_isValidJacobian)
				m_jacobianMatrix = g.m_jacobianMatrix;

			if (m_isValidGradient)
				m_gradientVector = g.m_gradientVector;

			if (m_isValidNewton)
				m_newtonVector = g.m_newtonVector;
		}

		return *this;
	}

	NOX::Abstract::Group&
	operator= (const NOX::Abstract::Group& g)
	{
		return operator=(dynamic_cast<const NOXGroup<T_Q,T_SIZE_MULTIVECT>&>(g));
	}

	void
	setX (const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& x)
	{
		resetIsValid();
		m_xVector = x;
	}

	void
	setX (const NOX::Abstract::Vector& x)
	{
		setX(dynamic_cast<const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>&>(x));
	}

	void
	computeX (const NOXGroup<T_Q,T_SIZE_MULTIVECT>& grp, const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& d, double step)
	{
		resetIsValid();
		m_xVector = grp.m_xVector + step*d;
	}

	void
	computeX (const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step)
	{
		const NOXGroup<T_Q,T_SIZE_MULTIVECT>& g = dynamic_cast<const NOXGroup<T_Q,T_SIZE_MULTIVECT>&>(grp);
		const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& v = dynamic_cast<const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>&>(d);
		computeX(g,v,step);
	}

	NOX::Abstract::Group::ReturnType
	computeF ()
	{
		if (m_isValidF)
			return NOX::Abstract::Group::Ok;

		m_isValidF = m_step.computeF(m_fVector,m_xVector);

		return (m_isValidF) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
	}

	NOX::Abstract::Group::ReturnType
	computeJacobian ()
	{
		if (m_isValidJacobian)
			return NOX::Abstract::Group::Ok;

		m_isValidJacobian = m_step.computeJacobian(m_jacobianMatrix,m_xVector);

		return (m_isValidJacobian) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
	}

	NOX::Abstract::Group::ReturnType
	computeGradient ()
	{
		if (m_isValidGradient)
			return NOX::Abstract::Group::Ok;

		if (!m_isValidF) {
			std::cerr << "ERROR: NOXGroup::computeGradient() - F is out of date wrt X" << std::endl;
			return NOX::Abstract::Group::BadDependency;
		}

		if (!m_isValidJacobian) {
			std::cerr << "ERROR: NOXGroup::computeGradient() - Jacobian is out of date wrt X" << std::endl;
			return NOX::Abstract::Group::BadDependency;
		}

		m_gradientVector = m_jacobianMatrix.transpose()*m_fVector;
		m_isValidGradient = true;

		return NOX::Abstract::Group::Ok;
	}

	NOX::Abstract::Group::ReturnType
	computeNewton (Teuchos::ParameterList& p)
	{
		if (m_isValidNewton)
			return NOX::Abstract::Group::Ok;

		if (!m_isValidF) {
			std::cerr << "ERROR: NOXGroup::computeNewton() - invalid F" << std::endl;
			throw "NOX Error";
		}

		if (!m_isValidJacobian) {
			std::cerr << "ERROR: NOXGroup::computeNewton() - invalid Jacobian" << std::endl;
			throw "NOX Error";
		}

		NOX::Abstract::Group::ReturnType status = applyJacobianInverse(p,m_fVector,m_newtonVector);
		m_isValidNewton = (status == NOX::Abstract::Group::Ok);

		m_newtonVector = -m_newtonVector;

		return status;
	}

	NOX::Abstract::Group::ReturnType
	applyJacobian (const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& input, NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& result) const
	{
		if (!m_isValidJacobian)
			return NOX::Abstract::Group::BadDependency;

		result = m_jacobianMatrix*input;

		return NOX::Abstract::Group::Ok;
	}

	NOX::Abstract::Group::ReturnType
	applyJacobian (const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
	{
		const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& v = dynamic_cast<const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>&>(input);
		NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& res = dynamic_cast<NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>&>(result);
		return applyJacobian(v,res);
	}

	NOX::Abstract::Group::ReturnType
	applyJacobianTranspose (const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& input, NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& result) const
	{
		if (!m_isValidJacobian)
			return NOX::Abstract::Group::BadDependency;

		result = m_jacobianMatrix.transpose()*input;

		return NOX::Abstract::Group::Ok;
	}

	NOX::Abstract::Group::ReturnType
	applyJacobianTranspose (const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
	{
		const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& v = dynamic_cast<const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>&>(input);
		NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& res = dynamic_cast<NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>&>(result);
		return applyJacobianTranspose(v,res);
	}

	NOX::Abstract::Group::ReturnType
	applyJacobianInverse (Teuchos::ParameterList& p, const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& input, NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& result) const
	{
		if (!m_isValidJacobian) {
			std::cerr << "ERROR: NOXGroup::applyJacobianInverse() - invalid Jacobian" << std::endl;
			throw "NOX Error";
		}

		// TODO : tester si le jacobien est inversible
		result = m_jacobianMatrix.inverse()*input;
		
		return NOX::Abstract::Group::Ok;
	}

	NOX::Abstract::Group::ReturnType
	applyJacobianInverse (Teuchos::ParameterList& p, const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const
	{
		const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& v = dynamic_cast<const NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>&>(input);
		NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>& res = dynamic_cast<NOXVector<T_Q::DOF*T_SIZE_MULTIVECT>&>(result);
		return applyJacobianInverse(p,v,res);
	}

	// TODO : Overloader applyJacobianInverseMultiVector avec une meilleure implementation que celle par defaut
	// (qui appelle plusieurs fois applyJacobianInverse)
	
	bool
	isF () const
	{
		return m_isValidF;
	}
	
	bool
	isJacobian () const
	{
		return m_isValidJacobian;
	}
	
	bool
	isGradient () const
	{
		return m_isValidGradient;
	}
	
	bool
	isNewton () const
	{
		return m_isValidNewton;
	}

	const NOX::Abstract::Vector&
	getX () const
	{
		return m_xVector;
	}

	const NOX::Abstract::Vector&
	getF () const
	{
		return m_fVector;
	}

	double
	getNormF () const
	{
		if (!m_isValidF) {
			std::cerr << "ERROR: NOXGroup::getNormF() - invalid F, please call computeF() first" << std::endl;
			throw "NOX Error";
		}

		return m_fVector.norm();
	}

	const NOX::Abstract::Vector&
	getGradient () const
	{
		return m_gradientVector;
	}

	const NOX::Abstract::Vector&
	getNewton () const
	{
		return m_newtonVector;
	}

	Teuchos::RCP<const NOX::Abstract::Vector>
	getXPtr () const
	{
		return Teuchos::RCP<const NOX::Abstract::Vector>(&m_xVector,false);
	}

	Teuchos::RCP<const NOX::Abstract::Vector>
	getFPtr () const
	{
		return Teuchos::RCP<const NOX::Abstract::Vector>(&m_fVector,false);
	}

	Teuchos::RCP<const NOX::Abstract::Vector>
	getGradientPtr () const
	{
		return Teuchos::RCP<const NOX::Abstract::Vector>(&m_gradientVector,false);
	}

	Teuchos::RCP<const NOX::Abstract::Vector>
	getNewtonPtr () const
	{
		return Teuchos::RCP<const NOX::Abstract::Vector>(&m_newtonVector,false);
	}

	Teuchos::RCP<NOX::Abstract::Group>
	clone (NOX::CopyType type = NOX::DeepCopy) const
	{
		Teuchos::RCP<NOX::Abstract::Group> tmp;
		tmp = Teuchos::rcp(new NOXGroup<T_Q,T_SIZE_MULTIVECT>(*this));
		return tmp;
	}

	void
	print () const
	{
		std::cout << "x = " << m_xVector << std::endl;

		if (m_isValidF) {
			std::cout << "F(x) = " << m_fVector << std::endl;
		} else {
			std::cout << "F(x) has not been computed" << std::endl;
		}
			
		std::cout << std::endl;
	}
protected:
	void
	resetIsValid ()
	{
		m_isValidF = false;
		m_isValidJacobian = false;
		m_isValidGradient = false;
		m_isValidNewton = false;
	}

protected:
	NOXVector<T_Q::DOF*T_SIZE_MULTIVECT> m_xVector;
	NOXVector<T_Q::DOF*T_SIZE_MULTIVECT> m_fVector;
	NOXVector<T_Q::DOF*T_SIZE_MULTIVECT> m_newtonVector;
	NOXVector<T_Q::DOF*T_SIZE_MULTIVECT> m_gradientVector;
	Eigen::Matrix<double,T_Q::DOF*T_SIZE_MULTIVECT,T_Q::DOF*T_SIZE_MULTIVECT> m_jacobianMatrix;

	Abstract::Step<T_Q,T_SIZE_MULTIVECT>& m_step;

	bool m_isValidF;
	bool m_isValidJacobian;
	bool m_isValidGradient;
	bool m_isValidNewton;
};

#endif
