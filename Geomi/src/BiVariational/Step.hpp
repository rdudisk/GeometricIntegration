#ifndef DEF_BIVARIATIONAL_STEP
#define DEF_BIVARIATIONAL_STEP

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

namespace BiVariational {

template <typename T_SCALAR,
		  typename T_Q,
		  typename T_VEL>
class Step
{
private:
	Teuchos::RCP<Teuchos::ParameterList>			m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests;
	Teuchos::RCP<NOX::Epetra::Group>				m_grp;
	Teuchos::RCP<NOX::Solver::Generic>				m_solver;
	Teuchos::RCP<Epetra_Vector>						m_soln;

	Abstract::LieProblem<T_SCALAR,T_Q,T_VEL>* m_problem;
	StepInternals<T_SCALAR,T_Q,T_VEL>* m_internals;
	Epetra_Comm* m_comm;

public:
	Step (Abstract::LieProblem<T_SCALAR,T_Q,T_VEL>& problem, Epetra_Comm& comm)
	:	m_problem(&problem),
		m_comm(&comm)
	{
		m_internals = new StepInternals<T_SCALAR,T_Q,T_VEL>(problem,comm);

		// NOX
		m_soln = m_internals->getInitialGuess();

		m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
		Teuchos::ParameterList& solverParameters = *(m_solverParametersPtr.get());
		solverParameters.set("Nonlinear Solver","Line Search Based");

		Teuchos::ParameterList& printParams = solverParameters.sublist("Printing");
		printParams.set("MyPID", m_comm->MyPID());
		printParams.set("Output Precision", 3);
		printParams.set("Output Processor", 0);
		bool verbose = true;
		if (verbose)
			printParams.set("Output Information",
					NOX::Utils::OuterIteration +
					NOX::Utils::OuterIterationStatusTest +
					NOX::Utils::InnerIteration +
					NOX::Utils::LinearSolverDetails +
					NOX::Utils::Parameters +
					NOX::Utils::Details +
					NOX::Utils::Warning +
					NOX::Utils::Debug +
					NOX::Utils::TestDetails +
					NOX::Utils::Error);
		else
			printParams.set("Output Information",
					NOX::Utils::TestDetails +
					NOX::Utils::Error);

		Teuchos::ParameterList& searchParameters = solverParameters.sublist("Line Search");
		searchParameters.set("Method","Full Step");

		Teuchos::ParameterList& directionParameters = solverParameters.sublist("Direction");
		directionParameters.set("Method","Newton");

		Teuchos::ParameterList& newtonParameters
			= directionParameters.sublist("Newton");
		newtonParameters.set("Forcing Term Method","Constant");

		Teuchos::ParameterList& linearSolverParameters =
			newtonParameters.sublist("Linear Solver");
		linearSolverParameters.set("Aztec Solver","GMRES");
		linearSolverParameters.set("Max Iterations",800);
		linearSolverParameters.set("Tolerance",1e-4);
		//linearSolverParameters.set("Preconditioner","New Ifpack");
		linearSolverParameters.set("Preconditioner","None");
		linearSolverParameters.set("Preconditioner Reuse Policy","Reuse");
		linearSolverParameters.set("Max Age Of Prec",5);

		// TODO: Error on this line
		//solverParameters.sublist("Solver Options").set("Status Test Check Type", NOX::StatusTest::Complete);

		Teuchos::RCP<StepInternals<T_SCALAR,T_Q,T_VEL>> internals_rcp = Teuchos::rcp(m_internals);

		// WARNING: m_soln était dans le code de la doc NOX "noxSol", un NOX::Epetra::Vector et non pas un Epetra_Vector
		Teuchos::RCP<NOX::Epetra::MatrixFree> MF
			= Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams,internals_rcp,*m_soln));

		Teuchos::RCP<NOX::Epetra::FiniteDifference> FD
			= Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams,internals_rcp,*m_soln));

		
		Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = internals_rcp;
		Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = internals_rcp; //MF;
		Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = FD;
		Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys
			= Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(
						printParams,
						linearSolverParameters,
						iJac, internals_rcp->jacobian(), iPrec, FD, *m_soln));
		/*
			= Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(
						printParams,
						linearSolverParameters,
						iJac, MF, iPrec, FD, *m_soln));
						*/

		NOX::Epetra::Vector initialGuess(m_soln,NOX::Epetra::Vector::CreateView);
		//Teuchos::RCP<NOX::Epetra::Group> grpPtr =
		this->m_grp = Teuchos::rcp(new NOX::Epetra::Group(printParams,iReq,initialGuess,linSys));

		/*
		Teuchos::RCP<NOX::StatusTest::NormF> absresid = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
		Teuchos::RCP<NOX::StatusTest::NormF> relresid = Teuchos::rcp(new NOX::StatusTest::NormF(*m_grp.get(),1.0e-2));
		Teuchos::RCP<NOX::StatusTest::NormUpdate> update = Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
		Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2,1.0e-8));
		Teuchos::RCP<NOX::StatusTest::Combo> converged = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
		converged->addStatusTest(absresid);
		converged->addStatusTest(relresid);
		converged->addStatusTest(update);
		converged->addStatusTest(wrms);
		Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
		Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = Teuchos::rcp(new NOX::StatusTest::FiniteValue);
		Teuchos::RCP<NOX::StatusTest::Combo> combo = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
		combo->addStatusTest(converged);
		combo->addStatusTest(maxiters);
		combo->addStatusTest(fv);
		this->m_statusTests = combo;
		*/
		// Alternative test
		Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,NOX::StatusTest::NormF::Unscaled));
		Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
		this->m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));

		this->m_solver = NOX::Solver::buildSolver(this->m_grp,this->m_statusTests,this->m_solverParametersPtr);
	}

	~Step ()
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
