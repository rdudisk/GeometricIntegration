#ifndef DEF_VARIATIONAL_STEP
#define DEF_VARIATIONAL_STEP

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

namespace Variational {

template <typename T_M,
		  typename T_Q,
		  typename T_TQ>
class Step : public Abstract::Step<T_M,T_Q,T_TQ>
{
private:
	Teuchos::RCP<Teuchos::ParameterList>			m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests;
	Teuchos::RCP<NOXGroup<T_Q,1>>					m_grp;
	Teuchos::RCP<NOX::Solver::Generic>				m_solver;

public:
	Step<T_M,T_Q,T_TQ> (Abstract::Problem<T_M,T_Q>& problem)
	:	Abstract::Step<T_M,T_Q,T_TQ>(problem)
	{
		// m_internals is already initalized by base constructor, overriding
		Abstract::StepInternals<T_M,T_Q,T_TQ>* internals
			= new Abstract::StepInternals<T_M,T_Q,T_TQ>(problem);
		this->m_internals = internals;

		// NOX
		this->m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
		Teuchos::ParameterList& solverParameters = *m_solverParametersPtr;
		solverParameters.set("Nonlinear Solver","Tensor Based");

		Teuchos::ParameterList& directionParameters
			= solverParameters.sublist("Direction");
		directionParameters.set("Method","Tensor");

		Teuchos::ParameterList& globalStrategyParameters
			= solverParameters.sublist("Line Search");
		globalStrategyParameters.set("Method","Curvilinear");

		Teuchos::ParameterList& lineSearchParameters
			= globalStrategyParameters.sublist(globalStrategyParameters.get("Method","Curvilinear"));
		lineSearchParameters.set("Lambda Selection","Halving");
		lineSearchParameters.set("Max Iters",20);

		Teuchos::RCP<NOX::StatusTest::NormF> statusTestA
			= Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,NOX::StatusTest::NormF::Unscaled));
		Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB
			= Teuchos::rcp(new NOX::StatusTest::MaxIters(50));
		this->m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));


		Abstract::StepInternals<T_M,T_Q,T_TQ>* tmp
			= static_cast<Abstract::StepInternals<T_M,T_Q,T_TQ>*>(this->m_internals);
		this->m_grp = Teuchos::rcp(new NOXGroup<T_Q,1>(*tmp));
		this->m_solver = NOX::Solver::buildSolver(this->m_grp,this->m_statusTests,this->m_solverParametersPtr);
	}

	~Step<T_M,T_Q,T_TQ> ()
	{ }

	const bool
	makeStep (void)
	{
		// TODO: a travailler !!
		T_Q* q1 = new T_Q();
		bool success = true;
		bool verbose = false;

		try {
			NOX::StatusTest::StatusType status = m_solver->solve();
			const NOXGroup<T_Q,1>& solnGrp = dynamic_cast<const NOXGroup<T_Q,1>&>(m_solver->getSolutionGroup());
			const NOXVector<T_Q::DOF>& solnVec = dynamic_cast<const NOXVector<T_Q::DOF>&>(solnGrp.getX());

			if (status != NOX::StatusTest::Converged)
				success = false;
		} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

		success &= this->m_internals->computeSolution();

		return success;
	}
};
} // namespace Variational

#endif
