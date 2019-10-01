#ifndef DEF_RKMK_DIAGONALSTEP
#define DEF_RKMK_DIAGONALSTEP

#include <NOX.H>
#include <NOX_Common.H>
#include <NOX_Utils.H>
#include <Teuchos_StandardCatchMacros.hpp>

namespace RKMK {

/**
 * This class inherits RKMK::Abstract::Step and implements a \f$p\f$ stages diagonal implicit RKMK step,
 * where \f$p\f$ is the template parameter T_N_INTERNAL_STAGES.
 * A step is diagonal implicit iff the coefficients of the Butcher tableau check \f$ a_{i,j}=0,\,\forall i\leq j\f$.
 *
 * The root \f$k_i\f$ of the equation \f$k_i-f\left(h\sum_{j=1}^ia_{i,j}k_j\right)=0\f$ is computed for every
 * ascending \f$i\f$ from \f$1\f$ to \f$p\f$ using the NOX library.
 * The configuration of the coefficients allows to solve for \f$p\f$ sets of \f$q\f$ equations
 * instead of one set of \f$p\times q\f$ equations, thus significantly decreasing the cost of the solver.
 */

template <typename T_LIE_ALGEBRA,
		  int T_N_INTERNAL_STAGES,
		  typename T_M = double>
class DiagonalStep : public Abstract::Step<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>
{
private:
	Teuchos::RCP<Teuchos::ParameterList>			m_solverParametersPtr;
	Teuchos::RCP<NOX::StatusTest::Combo>			m_statusTests;
	Teuchos::RCP<NOXGroup<NOXVector<T_LIE_ALGEBRA::DOF>,T_N_INTERNAL_STAGES>>
													m_grp;
	Teuchos::RCP<NOX::Solver::Generic>				m_solver;

public:
	DiagonalStep<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M> (Abstract::Problem<T_LIE_ALGEBRA,T_M>& problem)
	:	Abstract::Step<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>(problem)
	{
		// m_internals is already initalized by base constructor, overriding
		DiagonalStepInternals<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>* internals
			= new DiagonalStepInternals<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>(problem);
		this->m_internals = internals;

		// NOX
		this->m_solverParametersPtr = Teuchos::rcp(new Teuchos::ParameterList);
		Teuchos::ParameterList& solverParameters = *m_solverParametersPtr;
		solverParameters.set("Nonlinear Solver","Tensor Based");

		Teuchos::ParameterList& directionParameters = solverParameters.sublist("Direction");
		directionParameters.set("Method","Tensor");

		Teuchos::ParameterList& globalStrategyParameters = solverParameters.sublist("Line Search");
		globalStrategyParameters.set("Method","Curvilinear");

		Teuchos::ParameterList& lineSearchParameters = globalStrategyParameters.sublist(globalStrategyParameters.get("Method","Curvilinear"));
		lineSearchParameters.set("Lambda Selection","Halving");
		lineSearchParameters.set("Max Iters",20);

		Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,NOX::StatusTest::NormF::Unscaled));
		Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = Teuchos::rcp(new NOX::StatusTest::MaxIters(50));
		this->m_statusTests = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,statusTestA,statusTestB));


		DiagonalStepInternals<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>* tmp = static_cast<DiagonalStepInternals<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>*>(this->m_internals);
		this->m_grp = Teuchos::rcp(new NOXGroup<NOXVector<T_LIE_ALGEBRA::DOF>,T_N_INTERNAL_STAGES>(*tmp));
		this->m_solver = NOX::Solver::buildSolver(this->m_grp,this->m_statusTests,this->m_solverParametersPtr);
	}

	~DiagonalStep<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M> ()
	{ }

	const NOXVector<T_LIE_ALGEBRA::DOF>
	makeStep (void)
	{
		// TODO: faire des tests pour voir si tout se passe bien
		NOXVector<T_LIE_ALGEBRA::DOF>* Y1 = new NOXVector<T_LIE_ALGEBRA::DOF>();
		bool success = true;
		bool verbose = false;

		for (int i=0; i<T_N_INTERNAL_STAGES; i++) {
			try {
				(static_cast<DiagonalStepInternals<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>*>(this->m_internals))->currentStep(i);
				NOX::StatusTest::StatusType status = m_solver->solve();
				const NOXGroup<NOXVector<T_LIE_ALGEBRA::DOF>,1>& solnGrp = dynamic_cast<const NOXGroup<NOXVector<T_LIE_ALGEBRA::DOF>,1>&>(m_solver->getSolutionGroup());
				const NOXVector<T_LIE_ALGEBRA::DOF>& solnVec = dynamic_cast<const NOXVector<T_LIE_ALGEBRA::DOF>&>(solnGrp.getX());
				(static_cast<DiagonalStepInternals<T_LIE_ALGEBRA,T_N_INTERNAL_STAGES,T_M>*>(this->m_internals))->setSolutionK(solnVec);

				if (status != NOX::StatusTest::Converged)
					success = false;
			} TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
		}

		success &= this->m_internals->computeSolution();
		*Y1 = this->m_internals->reconstruct();

		return *Y1;
	}
};
} // namespace RKMK

#endif
