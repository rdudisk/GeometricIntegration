#ifndef DEF_VARIATIONAL_GALERKINSTEPINTERNALS
#define DEF_VARIATIONAL_GALERKINSTEPINTERNALS

#include <vector>

namespace Variational {

/**
 * Cette implémentation est un peu touffue, alors tâchons d'être méthodique.
 *
 * Le nombre de points de contrôle est désigné par \f$S+1\f$ dans les formules mathématiques.
 * Dans le code, \f$S\f$ est donné par T_N_STEPS.
 * L'indice associé est \f$\nu\in\{0,\ldots,S\}\f$ ou \f$\mu\f$.
 * Les points de contrôle sont les configurations \f$q^\nu\in Q\f$ où l'espace \f$Q\f$ possède T_Q::DOF degrés de liberté dénoté dans les formules par \f$N\f$.
 * L'ensemble des configurations associés à un pas de temps donné sont résumés par la "super"-configuration \f$\bar q\f$ de T_Q::DOF*(T_N_STEPS+1) degrés de liberté.
 * Il est représenté en mémoire par un std::vector<T_Q> et \f$q_\nu\f$ est accédé par bar_q[nu].
 *
 * Les polynômes de Lagrange \f$\phi_\nu\f$ associés aux points de contrôle \f$q^\nu\f$ sont de degré \f$S\f$ et définis pour \f$\alpha\in[0,1]\f$ par
 * \f[
 *		\phi_\nu(\alpha)=\prod_{0\leq\mu\leq S,\;\mu\neq\nu}\frac{\alpha-\alpha_\mu}{\alpha_\nu-\alpha_\mu}.
 * \f]
 * La courbe interpolée \f$q_d\f$ et sa dérivée \f$\dot q_d\f$ sont données en \f$t\in[0,h]\f$ respectivement par
 * \f[
 *		q_d(t;\bar q) = \sum_{0\leq\nu\leq S} q^\nu\phi_\nu(t/h),\qquad \dot q_d(t;\bar q) = \frac{1}{h}\sum_{0\leq\nu\leq S} q^\nu\dot\phi_\nu(t/h).
 * \f]
 *
 * Soit une quadrature \f$(c_k,w_k)\f$ d'ordre \f$r\f$ parcourue par l'indice \f$k\in\{0,\ldots,r-1\}\f$.
 * L'action est donnée pour \f$\bar q\f$ par
 * \f[
 *		\mathcal A_d(\bar q) = h\sum_{0\leq k<r}w_k L(c_kh;\bar q)
 * \f]
 * où \f$L(c_kh;\bar q) = L(q_d(c_kh;\bar q),\dot q_d(c_kh;\bar q))\f$ est un raccourci de notation.
 *
 * On note \f$F\f$ la fonction que l'on cherche à annuler.
 * Elle prend en argument les \f$S\f$ configurations \f$q_1,\ldots,q_S\f$ et retourne \f$S\f$ systèmes de \f$N\f$ équations.
 * Le premier système d'équations correspond à la DEL et les \f$S-1\f$ suivants aux conditions internes d'extrémalisation.
 * On notera l'ensemble des configurations \f$q_1,\ldots,q_S\f$ par \f$\bar q^*\f$, c'est à dire \f$\bar q\f$ sans la configuration \f$ q_0\f$, et
 * \f[
 *		F(\bar q^*)=\begin{pmatrix}F_0(\bar q^*)\\ \vdots \\ F_{S-1}(\bar q^*)\end{pmatrix}
 * \f]
 * où chaque \f$F_i(\bar q^*)\f$ est un vecteur colonne de taille \f$N\f$.
 *
 * Soit \f$\bar q_{n-1}\f$ et \f$\bar q_n\f$ les ensembles de configurations de deux pas de temps consécutifs,
 * la DEL est donnée par
 * \f[
 *		F_0(\bar q_n^*) = \sum_{0\leq k<r}w_k\left(h\left(\phi_0(c_k)\frac{\partial L}{\partial q}(c_kh;\bar q_n)+\phi_s(c_k)\frac{\partial L}{\partial q}(c_kh;\bar q_{n-1})\right)
 *			+\left(\dot\phi_0(c_k)\frac{\partial L}{\partial\dot q}(c_kh;\bar q_n)+\dot \phi_s(c_k)\frac{\partial L}{\partial\dot q}(c_kh;\bar q_{n-1})\right)\right) = 0
 * \f]
 * et les conditions internes d'extrémalisation d'action sont données pour tout \f$\nu\in\{1,\ldots,S-1\}\f$ par
 * \f[
 *		F_\nu(\bar q_n^*)=\sum_{0\leq k<r}w_k\left(h\phi_\nu(c_k)\frac{\partial L}{\partial q}(c_kh;\bar q_n)+\dot\phi_\nu(c_k)\frac{\partial L}{\partial \dot q}(c_kh;\bar q_n)\right)=0.
 * \f]
 *
 * Le jacobien de \f$F\f$ est donné par
 * \f[
 *		J(\bar q^*)=\left(\frac{\partial F_i}{\partial q_{j+1}}(\bar q^*)\right)_{ij}.
 * \f]
 * Pour \f$0\leq i\leq S-1\f$, \f$1\leq j\leq S\f$, et \f$\bar q=\bar q_n\f$, on calcule
 * \f[
 *		\frac{\partial F_i}{\partial q_j}(\bar q^*) = \sum_{0\leq k<r}w_k\left(\phi_i(c_k)\left(h\phi_j(c_k)\frac{\partial}{\partial q}\frac{\partial L}{\partial q}
 *																									+\dot\phi_j(c_k)\frac{\partial}{\partial v}\frac{\partial L}{\partial q}\right)
 *																				+\dot\phi_i(c_k)\left(\phi_j(c_k)\frac{\partial}{\partial q}\frac{\partial L}{\partial v}
 *																									+\frac{1}{h}\dot\phi_j(c_k)\frac{\partial}{\partial v}\frac{\partial L}{\partial v}\right)\right)
 * \f]
 * où toutes les dérivées secondes du lagrangien sont prises en \f$(c_k;\bar q_n)\f$.
 * On vérifie en particulier que la différentielle de DEL donne bien
 * \f[
 *		\frac{\partial F_0}{\partial q_j}(\bar q^*) = \sum_{0\leq k<r}w_k\left(\phi_0(c_k)\left(h\phi_j(c_k)\frac{\partial}{\partial q}\frac{\partial L}{\partial q}
 *																									+\dot\phi_j(c_k)\frac{\partial}{\partial v}\frac{\partial L}{\partial q}\right)
 *																				+\dot\phi_0(c_k)\left(\phi_j(c_k)\frac{\partial}{\partial q}\frac{\partial L}{\partial v}
 *																									+\frac{1}{h}\dot\phi_j(c_k)\frac{\partial}{\partial v}\frac{\partial L}{\partial v}\right)\right)
 * \f]
 */
template <typename T_M,
		  typename T_Q,
		  typename T_TQ,
		  int T_N_STEPS>
class GalerkinStepInternals : public Abstract::StepInternals<T_M,T_Q,T_TQ>, public ::Abstract::NOXStep<T_Q,T_N_STEPS>
{
protected:
	/** \f$\bar q_n\f$, de longueur T_N_STEPS+1 */
	std::vector<T_Q>		m_v_cur_q;
	/** \f$\bar q_{n-1}\f$, de longueur T_N_STEPS+1 */
	std::vector<T_Q>		m_v_prev_q;
	LagrangeInterpolation<T_Q>	m_interp;
	int						m_quad_deg;

public:
	GalerkinStepInternals<T_M,T_Q,T_TQ,T_N_STEPS> (Abstract::Problem<T_M,T_Q>& problem, int quad_deg)
	:	Abstract::StepInternals<T_M,T_Q,T_TQ>(problem)
	{
		m_interp = LagrangeInterpolation<T_Q>();
		m_quad_deg = quad_deg;
	}

	/*
	 * T_Q
	posFromVel (T_Q q0, T_TQ v0)
	{
		return q0+this->m_h*v0;
	}
	*/

	// TODO
	const NOXVector<T_Q::DOF*T_N_STEPS>&
	getInitialGuess ()
	{
		NOXVector<T_Q::DOF*T_N_STEPS>* ret;
		T_Q q0 = m_v_prev_q[0];
		T_Q q1 = m_v_prev_q[T_N_STEPS];

		for (int i=1; i<=T_N_STEPS; i++) {
			ret->segment(T_Q::DOF,T_N_STEPS) = NOXVector<T_Q::DOF>(q1+(float(i)/(T_N_STEPS*this->m_h))*(q1-q0));
		}
		return *ret;
	}

	void
	updatePosition (const NOXVector<T_Q::DOF*T_N_STEPS>& q)
	{
		T_Q tmp = m_v_prev_q[T_N_STEPS];
		m_v_prev_q.clear();
		m_v_prev_q.push_back(tmp);

		for (int i=1; i<=T_N_STEPS; i++) {
			m_v_prev_q.push_back(T_Q(q.segment(T_Q::DOF,i*T_Q::DOF)));
		}
	}

	bool
	computeF (NOXVector<T_Q::DOF*T_N_STEPS>& f, const NOXVector<T_Q::DOF*T_N_STEPS>& q)
	{
		int nu;	// index polynome
		int k;	// index date


		// update current position set
		m_v_cur_q.clear();
		m_v_cur_q.push_back(m_v_prev_q[T_N_STEPS]);
		for (nu=0; nu<T_N_STEPS; nu++) {
			m_v_cur_q.push_back(T_Q(q.segment(T_Q::DOF,nu*T_Q::DOF)));
		}

		int r = m_quad_deg;		// degré quadrature
		std::vector<double> w = GaussLegendre::weights(r);	// poids quadrature, len = r
		std::vector<double> c = GaussLegendre::dates(r);		// dates quadrature, len = r


		//  pour chaque date d'indice k, contient le vecteur des T_N_STEPS polynomes de lagrange pris en ce point
		std::vector<std::vector<double>> vv_lag;
		// idem pour les dérivées
		std::vector<std::vector<double>> vv_lag_der;

		vv_lag = this->m_interp.polynomials(T_N_STEPS,c);
		vv_lag_der = this->m_interp.polynomials_derivatives(T_N_STEPS,c);


		// vecteur des r positions correspondant aux dates c pour la paramétrisation m_qCurrent;
		std::vector<T_Q> v_pos_interp = m_interp.pos_interp(T_N_STEPS,c,m_v_cur_q);
		// idem pour les vitesses
		std::vector<T_Q> v_vel_interp = m_interp.vel_interp(T_N_STEPS,c,m_v_cur_q,this->m_h);
		// vecteur des r positions correspondant aux dates c pour la paramétrisation m_qPrevious;
		std::vector<T_Q> v_prev_pos_interp = m_interp.pos_interp(T_N_STEPS,c,m_v_prev_q);
		// idem pour les vitesses
		std::vector<T_Q> v_prev_vel_interp = m_interp.vel_interp(T_N_STEPS,c,m_v_prev_q,this->m_h);
		
		// vecteur des r différentielles du lagrangien par rapport à q prises aux r dates c
		std::vector<NOXVector<T_Q::DOF>> v_cur_dLdq;
		// idem pour les différentielles par rapport à v
		std::vector<NOXVector<T_Q::DOF>> v_cur_dLdv;
		// idem pour le pas de temps précédent 
		std::vector<NOXVector<T_Q::DOF>> v_prev_dLdq;
		// idem par rapport à v
		std::vector<NOXVector<T_Q::DOF>> v_prev_dLdv;

		for (k=0; k<r; k++) {
			v_cur_dLdq.push_back(this->m_problem.dLdq(v_pos_interp[k],v_vel_interp[k]));
			v_cur_dLdv.push_back(this->m_problem.dLdv(v_pos_interp[k],v_vel_interp[k]));
			v_prev_dLdq.push_back(this->m_problem.dLdq(v_prev_pos_interp[k],v_prev_vel_interp[k]));
			v_prev_dLdv.push_back(this->m_problem.dLdv(v_prev_pos_interp[k],v_prev_vel_interp[k]));
		}

		f = NOXVector<T_Q::DOF*T_N_STEPS>::Zero();

		// DEL
		for (k=0; k<r; k++) {
			f.head(T_Q::DOF) += w[k]*(this->m_h*vv_lag[k][0]*v_cur_dLdq[k]+vv_lag_der[k][0]*v_cur_dLdv[k]+this->m_h*vv_lag[k][T_N_STEPS]*v_prev_dLdq[k]+vv_lag_der[k][T_N_STEPS]*v_prev_dLdv[k]);
		}

		NOXVector<T_Q::DOF> somme;
		// Internal equations
		for (nu=1;nu<T_N_STEPS;nu++) {
			somme = T_Q::Zero();
			for (k=0;k<r;k++) {
				// les indices sont faux, corriger
				somme += w[k]*(this->m_h*vv_lag[k][nu]*v_cur_dLdq[k]+vv_lag_der[k][nu]*v_cur_dLdv[k]);
			}
			f.segment(T_Q::DOF,nu*T_Q::DOF) = somme;
		}

		return true;
	}

	bool
	computeInitF (NOXVector<T_Q::DOF*(T_N_STEPS-1)>& f, const NOXVector<T_Q::DOF*(T_N_STEPS-1)>& q)
	{
		int nu;	// index polynome
		int k;	// index date


		// update current position set
		// don't touch q[0] and q[T_N_STEPS]
		for (nu=0; nu<T_N_STEPS-1; nu++) {
			m_v_cur_q[nu+1] = T_Q(q.segment(T_Q::DOF,nu*T_Q::DOF));
		}

		int r = m_quad_deg;		// degré quadrature
		std::vector<double> w = GaussLegendre::weights(r);	// poids quadrature, len = r
		std::vector<double> c = GaussLegendre::dates(r);		// dates quadrature, len = r


		//  pour chaque date d'indice k, contient le vecteur des T_N_STEPS polynomes de lagrange pris en ce point
		std::vector<std::vector<double>> vv_lag;
		// idem pour les dérivées
		std::vector<std::vector<double>> vv_lag_der;

		vv_lag = this->m_interp.polynomials(T_N_STEPS,c);
		vv_lag_der = this->m_interp.polynomials_derivatives(T_N_STEPS,c);


		// vecteur des r positions correspondant aux dates c pour la paramétrisation m_qCurrent;
		std::vector<T_Q> v_pos_interp = m_interp.pos_interp(T_N_STEPS,c,m_v_cur_q);
		// idem pour les vitesses
		std::vector<T_Q> v_vel_interp = m_interp.vel_interp(T_N_STEPS,c,m_v_cur_q,this->m_h);
		
		NOXVector<T_Q::DOF> somme();
		// vecteur des r différentielles du lagrangien par rapport à q prises aux r dates c
		std::vector<NOXVector<T_Q::DOF>> v_cur_dLdq;
		// idem pour les différentielles par rapport à v
		std::vector<NOXVector<T_Q::DOF>> v_cur_dLdv;

		for (k=0; k<r; k++) {
			v_cur_dLdq.push_back(this->m_problem.dLdq(v_pos_interp[k],v_vel_interp[k]));
			v_cur_dLdv.push_back(this->m_problem.dLdv(v_pos_interp[k],v_vel_interp[k]));
		}

		// Internal equations
		for (nu=1;nu<T_N_STEPS;nu++) {
			somme = T_Q::Zero();
			for (k=0;k<r;k++) {
				// les indices sont faux, corriger
				somme += w[k]*(this->m_h*vv_lag[k][nu]*v_cur_dLdq[k]+vv_lag_der[k][nu]*v_cur_dLdv[k]);
			}
			f.segment(T_Q::DOF,(nu-1)*T_Q::DOF) = somme;
		}

		return true;
	}

	/**
	 * \f[
	 *		J_{ij} = \frac{\partial F_i}{\partial q_j}
	 *	\f]
	 */
	bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF*T_N_STEPS,T_Q::DOF*T_N_STEPS>& J, const NOXVector<T_Q::DOF*T_N_STEPS>& q)
	{
		int nu;	// index polynome
		int k;	// index date
		int i;	// jacobian line index
		int j;	// jacobian column index


		// update current position set
		m_v_cur_q.clear();
		m_v_cur_q.push_back(m_v_prev_q[T_N_STEPS]);
		for (nu=0; nu<T_N_STEPS; nu++) {
			m_v_cur_q.push_back(T_Q(q.segment(T_Q::DOF,nu)));
		}


		int r = m_quad_deg;		// degré quadrature
		std::vector<double> w = GaussLegendre::weights(r);	// poids quadrature, len = r
		std::vector<double> c = GaussLegendre::dates(r);		// dates quadrature, len = r


		//  pour chaque date d'indice k, contient le vecteur des T_N_STEPS polynomes de lagrange pris en ce point
		std::vector<std::vector<double>> vv_lag;
		// idem pour les dérivées
		std::vector<std::vector<double>> vv_lag_der;

		vv_lag = this->m_interp.polynomials(T_N_STEPS,c);
		vv_lag_der = this->m_interp.polynomials_derivatives(T_N_STEPS,c);


		// vecteur des r positions correspondant aux dates c pour la paramétrisation m_qCurrent;
		std::vector<T_Q> v_pos_interp = m_interp.pos_interp(T_N_STEPS,c,m_v_cur_q);
		// idem pour les vitesses
		std::vector<T_Q> v_vel_interp = m_interp.vel_interp(T_N_STEPS,c,m_v_cur_q,this->m_h);
		// vecteur des r positions correspondant aux dates c pour la paramétrisation m_qPrevious;
		std::vector<T_Q> v_prev_pos_interp = m_interp.pos_interp(T_N_STEPS,c,m_v_prev_q);
		// idem pour les vitesses
		std::vector<T_Q> v_prev_vel_interp = m_interp.vel_interp(T_N_STEPS,c,m_v_prev_q,this->m_h);
		
		// vecteur des r jacobiens du lagrangien par rapport à qq prises aux r dates c
		std::vector<Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>> v_JqdLdq;
		// idem pour qv
		std::vector<Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>> v_JqdLdv;
		// idem pour vq
		std::vector<Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>> v_JvdLdq;
		// idem pour vv
		std::vector<Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>> v_JvdLdv;

		for (k=0; k<r; k++) {
			v_JqdLdq.push_back(this->m_problem.JqdLdq(v_pos_interp[k],v_vel_interp[k]));
			v_JqdLdv.push_back(this->m_problem.JqdLdv(v_pos_interp[k],v_vel_interp[k]));
			v_JvdLdq.push_back(this->m_problem.JvdLdq(v_pos_interp[k],v_vel_interp[k]));
			v_JvdLdv.push_back(this->m_problem.JvdLdv(v_pos_interp[k],v_vel_interp[k]));
		}

		//J = Eigen::Matrix<double,T_Q::DOF*T_N_STEPS,T_Q::DOF*T_N_STEPS>::Zero();
		Eigen::Matrix<double,T_Q::DOF,T_Q::DOF> somme = Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>::Zero();

		for (j=0; j<T_N_STEPS; j++) {
			/*somme =  Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>::Zero();
			for (k=0; k<r; k++) {
				somme += w[k]*(vv_lag[k][0]*(this->m_h*vv_lag[k][j]*v_JqdLdq[k]+vv_lag_der[k][j]*v_JvdLdq[k]) + vv_lag_der[k][0]*(v_lag[k][j]*v_JqdLdv[k]+(vv_lag_der[k][j]/this->m_h)*v_JvdLdv[k]));
			}
			J.block<T_Q::DOF,T_Q::DOF>(0,j*T_Q::DOF) = somme;*/
			for (i=0; i<T_N_STEPS; i++) {
				somme =  Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>::Zero();
				for (k=0; k<r; k++) {
					somme += w[k]*(vv_lag[k][i]*(this->m_h*vv_lag[k][j]*v_JqdLdq[k]+vv_lag_der[k][j]*v_JvdLdq[k]) + vv_lag_der[k][i]*(vv_lag[k][j]*v_JqdLdv[k]+(vv_lag_der[k][j]/this->m_h)*v_JvdLdv[k]));
				}
				J.block(T_Q::DOF,T_Q::DOF,i*T_Q::DOF,j*T_Q::DOF) = somme;
			}
		}

		return true;
	}
	
	bool
	computeInitJacobian (Eigen::Matrix<double,T_Q::DOF*T_N_STEPS,T_Q::DOF*T_N_STEPS>& J, const NOXVector<T_Q::DOF*T_N_STEPS>& q)
	{
		int nu;	// index polynome
		int k;	// index date
		int i;	// jacobian line index
		int j;	// jacobian column index


		// update current position set
		// don't touch q[0] and q[T_N_STEPS]
		for (nu=0; nu<T_N_STEPS-1; nu++) {
			m_v_cur_q[nu+1] = T_Q(q.segment(T_Q::DOF,nu*T_Q::DOF));
		}


		int r = m_quad_deg;		// degré quadrature
		std::vector<double> w = GaussLegendre::weights(r);	// poids quadrature, len = r
		std::vector<double> c = GaussLegendre::dates(r);		// dates quadrature, len = r


		//  pour chaque date d'indice k, contient le vecteur des T_N_STEPS polynomes de lagrange pris en ce point
		std::vector<std::vector<double>> vv_lag;
		// idem pour les dérivées
		std::vector<std::vector<double>> vv_lag_der;

		vv_lag = this->m_interp.polynomials(T_N_STEPS,c);
		vv_lag_der = this->m_interp.polynomials_derivatives(T_N_STEPS,c);


		// vecteur des r positions correspondant aux dates c pour la paramétrisation m_qCurrent;
		std::vector<T_Q> v_pos_interp = m_interp.pos_interp(T_N_STEPS,c,m_v_cur_q);
		// idem pour les vitesses
		std::vector<T_Q> v_vel_interp = m_interp.vel_interp(T_N_STEPS,c,m_v_cur_q,this->m_h);
		
		// vecteur des r jacobiens du lagrangien par rapport à qq prises aux r dates c
		std::vector<Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>> v_JqdLdq;
		// idem pour qv
		std::vector<Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>> v_JqdLdv;
		std::vector<Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>> v_JvdLdq;
		std::vector<Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>> v_JvdLdv;

		for (k=0; k<r; k++) {
			v_JqdLdq.push_back(this->m_problem.JqdLdq(v_pos_interp[k],v_vel_interp[k]));
			v_JqdLdv.push_back(this->m_problem.JqdLdv(v_pos_interp[k],v_vel_interp[k]));
			v_JvdLdq.push_back(this->m_problem.JvdLdq(v_pos_interp[k],v_vel_interp[k]));
			v_JvdLdv.push_back(this->m_problem.JvdLdv(v_pos_interp[k],v_vel_interp[k]));
		}

		//J = Eigen::Matrix<double,T_Q::DOF*T_N_STEPS,T_Q::DOF*T_N_STEPS>::Zero();
		Eigen::Matrix<double,T_Q::DOF,T_Q::DOF> somme = Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>::Zero();

		for (j=1; j<T_N_STEPS; j++) {
			/*somme =  Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>::Zero();
			for (k=0; k<r; k++) {
				somme += w[k]*(vv_lag[k][0]*(this->m_h*vv_lag[k][j]*v_JqdLdq[k]+vv_lag_der[k][j]*v_JvdLdq[k]) + vv_lag_der[k][0]*(v_lag[k][j]*v_JqdLdv[k]+(vv_lag_der[k][j]/this->m_h)*v_JvdLdv[k]));
			}
			J.block<T_Q::DOF,T_Q::DOF>(0,j*T_Q::DOF) = somme;*/
			for (i=1; i<T_N_STEPS; i++) {
				somme =  Eigen::Matrix<double,T_Q::DOF,T_Q::DOF>::Zero();
				for (k=0; k<r; k++) {
					somme += w[k]*(vv_lag[k][i]*(this->m_h*vv_lag[k][j]*v_JqdLdq[k]+vv_lag_der[k][j]*v_JvdLdq[k]) + vv_lag_der[k][i]*(vv_lag[k][j]*v_JqdLdv[k]+(vv_lag_der[k][j]/this->m_h)*v_JvdLdv[k]));
				}
				J.block<T_Q::DOF,T_Q::DOF>((i-1)*T_Q::DOF,(j-1)*T_Q::DOF) = somme;
			}
		}

		return true;
	}
};

template <typename T_M,
		  typename T_Q,
		  typename T_TQ,
		  int T_N_STEPS>
class GalerkinStepInitWrapper : public ::Abstract::NOXStep<T_Q,T_N_STEPS-1>
{
private:
	GalerkinStepInternals<T_M,T_Q,T_TQ,T_N_STEPS>* m_internals;

public:
	GalerkinStepInitWrapper<T_M,T_Q,T_TQ,T_N_STEPS> (GalerkinStepInternals<T_M,T_Q,T_TQ,T_N_STEPS>* internals)
	{
		m_internals = internals;
	}

	const NOXVector<T_Q::DOF*(T_N_STEPS-1)>&
	getInitialGuess ()
	{
		// TODO
	}

	bool
	computeF (NOXVector<T_Q::DOF*(T_N_STEPS-1)>& f, const NOXVector<T_Q::DOF*(T_N_STEPS-1)>& q)
	{
		return m_internals->computeInitF(f,q);
	}

	bool
	computeJacobian (Eigen::Matrix<double,T_Q::DOF*(T_N_STEPS-1),T_Q::DOF*(T_N_STEPS-1)>& J, const NOXVector<T_Q::DOF*(T_N_STEPS-1)>& q)
	{
		return m_internals->computeInitJacobian(J,q);
	}
};

} // namespace Variational

#endif
