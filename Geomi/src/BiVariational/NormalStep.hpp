#ifndef DEF_BIVARIATIONAL_NORMAL_STEP
#define DEF_BIVARIATIONAL_NORMAL_STEP

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

namespace BiVariational {

template <typename T_SCALAR,
		  typename T_Q,
		  typename T_VEL>
class NormalStep
{
private:
	Teuchos::RCP<Teuchos::ParameterList>			m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests;
	Teuchos::RCP<NOX::Epetra::Group>				m_grp;
	Teuchos::RCP<NOX::Solver::Generic>				m_solver;
	//Teuchos::RCP<Epetra_Vector>						m_soln;

	Abstract::Problem<T_SCALAR,T_Q,T_VEL>* m_problem;
	NormalStepInternals<T_SCALAR,T_Q,T_VEL>* m_internals;

public:
	NormalStep (Abstract::Problem<T_SCALAR,T_Q,T_VEL>& problem)
	:	m_problem(&problem)
	{
		m_internals = new NormalStepInternals<T_SCALAR,T_Q,T_VEL>(problem);

		// NOX
		//m_soln = m_internals->getInitialGuess();

		m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
		Teuchos::ParameterList& solverParameters = *(m_solverParametersPtr.get());
		solverParameters.set("Nonlinear Solver","Line Search Based");

		Teuchos::ParameterList& directionParameters
			= solverParameters.sublist("Direction");
		directionParameters.set("Method","Tensor");

		Teuchos::ParameterList& globalStrategyParameters
			= solverParameters.sublist("Line Search");
		searchParameters.set("Method","Curvilinear");

		Teuchos::ParameterList& linearSearchParameters
			= globalStrategyParameters.sublist(globalStrategyParameters.get("Method","Curvilinear"));
		lineSearchParameters.set("Lambda Selection","Halving");
		lineSearchParameters.set("Max Iters",100);

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

	~NormalStep ()
	{
		delete m_internals;
	}

	bool
	initialize ( )
	{
		bool success = m_internals->initLimit();
		return success;
	}

	const bool
	makeStep (void)
	{
		bool success = true;
		bool verbose = false;

		try {
			m_internals->setupStep();
			m_solver->reset(NOX::Epetra::Vector(this->m_internals->getInitialGuess()));
			NOX::StatusTest::StatusType status = m_solver->solve();
			const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(m_solver->getSolutionGroup());
			const Epetra_Vector& solution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();


			//std::cout << "Approximate solution:" << std::endl << solution  << std::endl;

			if (status != NOX::StatusTest::Converged)
				success = false;
			
			m_internals->setSolution(solution);
		} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

		std::cout << "makeStep success: " << success << std::endl;

		//success &= this->m_internals->computeSolution();

		// TODO: après la résolution, mettre à jour l'indice temporel global
		// pour les internals

		return success;
	}
};
} // namespace Variational

#endif
