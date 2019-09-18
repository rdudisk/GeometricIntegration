#ifndef DEF_VARIATIONAL_COVARIANTSTEP
#define DEF_VARIATIONAL_COVARIANTSTEP

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

namespace Variational {

template <typename T_M,
		  typename T_Q,
		  typename T_ALGEBRA>
class CovariantStep : public Abstract::Step<T_M,
											T_Q,
											CovariantStepInternals<T_M,T_Q,T_ALGEBRA>,
											Abstract::LieProblem<T_M,T_Q,T_ALGEBRA>>
{
	using Internals = CovariantStepInternals<T_M,T_Q,T_ALGEBRA>;
	using Problem = Abstract::LieProblem<T_M,T_Q,T_ALGEBRA>;

private:
	Teuchos::RCP<Teuchos::ParameterList>			m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests;
	Teuchos::RCP<NOXGroup<T_ALGEBRA,1>>				m_grp;
	Teuchos::RCP<NOX::Solver::Generic>				m_solver;

public:
	CovariantStep<T_M,T_Q,T_ALGEBRA> (Problem& problem)
	:	Abstract::Step<T_M,T_Q,Internals,Problem>(problem)
	{
		// NOX
		this->m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
		Teuchos::ParameterList& solverParameters = *m_solverParametersPtr;

		solverParameters.set("Nonlinear Solver","Line Search Based");
		Teuchos::ParameterList& lineSearchParameters = solverParameters.sublist("Line Search");
		lineSearchParameters.set("Method","Full Step");

		Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,NOX::StatusTest::NormF::Unscaled));
		Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(500));
		this->m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));

		this->m_grp = Teuchos::rcp(new NOXGroup<T_Q,1>(*(this->m_internals)));
		this->m_solver = NOX::Solver::buildSolver(this->m_grp,this->m_statusTests,this->m_solverParametersPtr);
	}

	~CovariantStep<T_M,T_Q,T_ALGEBRA> ()
	{ }

	void
	setData (T_M h_var, T_Q q0_var, T_Q q1_var)
	{ this->m_internals->setData(h_var,q0_var,q1_var); }
	
	void
	initialize (void)
	{ }

	const T_Q
	makeStep (void)
	{
		bool success = true;
		bool verbose = false;

		try {
			m_solver->reset(this->m_internals->getInitialGuess());
			NOX::StatusTest::StatusType status = m_solver->solve();
			const NOXGroup<T_ALGEBRA,1>& solnGrp = dynamic_cast<const NOXGroup<T_ALGEBRA,1>&>(m_solver->getSolutionGroup());
			const NOXVector<T_ALGEBRA::DOF>& solnVec = dynamic_cast<const NOXVector<T_ALGEBRA::DOF>&>(solnGrp.getX());

			if (status != NOX::StatusTest::Converged)
				success = false;

			// TODO: la ligne suivante n'est pas bonne.
			// Il faut faire l'update par exponentielle
			T_Q solution(solnVec);

			return solution;
		} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
	}
};
} // namespace Variational

#endif
