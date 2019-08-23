#ifndef DEF_VARIATIONAL_MIDPOINTSTEP
#define DEF_VARIATIONAL_MIDPOINTSTEP

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

namespace Variational {

template <typename T_M,
		  typename T_Q,
		  typename T_TQ>
class MidpointStep : public Abstract::Step<T_M,T_Q,T_TQ>
{
private:
	Teuchos::RCP<Teuchos::ParameterList>			m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests;
	Teuchos::RCP<NOXGroup<T_Q,1>>					m_grp;
	Teuchos::RCP<NOX::Solver::Generic>				m_solver;

public:
	MidpointStep<T_M,T_Q,T_TQ> (Abstract::Problem<T_M,T_Q>& problem)
	:	Abstract::Step<T_M,T_Q,T_TQ>(problem)
	{
		// m_internals is already initalized by base constructor, overriding
		MidpointStepInternals<T_M,T_Q,T_TQ>* internals = new MidpointStepInternals<T_M,T_Q,T_TQ>(problem);
		this->m_internals = internals;

		// NOX
		this->m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
		Teuchos::ParameterList& solverParameters = *m_solverParametersPtr;

		solverParameters.set("Nonlinear Solver","Line Search Based");
		Teuchos::ParameterList& lineSearchParameters = solverParameters.sublist("Line Search");
		lineSearchParameters.set("Method","Full Step");

		Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,NOX::StatusTest::NormF::Unscaled));
		Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(500));
		this->m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));

		MidpointStepInternals<T_M,T_Q,T_TQ>* tmp = static_cast<MidpointStepInternals<T_M,T_Q,T_TQ>*>(this->m_internals);
		this->m_grp = Teuchos::rcp(new NOXGroup<T_Q,1>(*tmp));
		this->m_solver = NOX::Solver::buildSolver(this->m_grp,this->m_statusTests,this->m_solverParametersPtr);
	}

	~MidpointStep<T_M,T_Q,T_TQ> ()
	{ }

	const T_Q
	makeStep (void)
	{
		// TODO: faire des tests pour voir si tout se passe bien
		bool success = true;
		bool verbose = false;

		try {
			m_solver->reset(static_cast<MidpointStepInternals<T_M,T_Q,T_TQ>*>(this->m_internals)->getInitialGuess());
			NOX::StatusTest::StatusType status = m_solver->solve();
			const NOXGroup<T_Q,1>& solnGrp = dynamic_cast<const NOXGroup<T_Q,1>&>(m_solver->getSolutionGroup());
			const NOXVector<T_Q::DOF>& solnVec = dynamic_cast<const NOXVector<T_Q::DOF>&>(solnGrp.getX());

			if (status != NOX::StatusTest::Converged)
				success = false;

			T_Q solution(solnVec);

			return solution;
		} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
	}
};
} // namespace Variational

#endif
