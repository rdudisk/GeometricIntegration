#ifndef DEF_COMMON_LEGENDREINTERPOLATION
#define DEF_COMMON_LEGENDREINTERPOLATION

#include <vector>

class LagrangePolynomials
{
protected:
	std::vector<double> m_dates;
	int m_degree;
	bool m_computed;
	std::vector<std::vector<double>> m_polynomials;
	std::vector<std::vector<double>> m_polynomials_derivatives;

public:
	LagrangePolynomials ()
	{
		m_degree = 0;
		m_dates.clear();
		m_computed = false;
	}

	std::vector<double>
	chebyshev (int degree)
	{
		double den = 1.0/double(degree);
		std::vector<double> ret;
		for (int i=0;i<=degree;i++) {
			ret.push_back(double(i)*den);
		}
		return ret;
	}

	void
	check_computed (int s, std::vector<double> dates)
	{
		if (s!=m_degree) {
			m_degree = s;
			m_computed = false;
		}
		if (dates!=m_dates) {
			m_dates = dates;
			m_computed = false;
		}
		if (!m_computed) {
			compute_polynomials();
			compute_polynomials_derivatives();
			m_computed = true;
		}
	}

	std::vector<std::vector<double>>
	polynomials (int s, std::vector<double> dates)
	{
		check_computed(s, dates);
		return m_polynomials;
	}

	std::vector<std::vector<double>>
	polynomials_derivatives (int s, std::vector<double> dates)
	{
		check_computed(s, dates);
		return m_polynomials_derivatives;
	}

	// méthode brutale, peut-être un peu d'optimisation serait bienvenue...
	/**
	 * Calcul les polynomes de Lagrange
	 * \f[
	 *		\phi_j(d_i)=\prod_{0\leq k\leq s,\;k\neq j}\frac{d_i-\alpha_k}{\alpha_j-\alpha_k}
	 * \f]
	 * Le résultat est un vecteur de taille m_dates.size() (typiquement le degré d'interpolation)
	 * contenant pour chaque date un vecteur des valeurs des m_degree polynomes de Lagrange (typiquement s=m_degree=T_N_STEPS).
	 * Plus concrètement, vec[i][j]=\f$\phi_j(d_i)\f$.
	 */
	void
	compute_polynomials ()
	{
		int i; // index date
		int j; // index polynome
		int k; // index produit
		//std::vector<std::vector<double>> ret;
		std::vector<double> tmp;
		std::vector<double> coeffs = chebyshev(m_degree);
		double prod;
		double t;

		m_polynomials.clear();
		for (i=0; i<m_dates.size(); i++) {
			t = m_dates[i];
			tmp.clear();
			for (j=0; j<=m_degree; j++) {
				prod = 1.0;
				for (k=0; k<=m_degree; k++) {
					if (j!=k)
						prod *= (t-coeffs[k])/(coeffs[j]-coeffs[k]);
				}
				tmp.push_back(prod);
			}
			m_polynomials.push_back(tmp);
		}
	}

	// méthode brutale, peut-être un peu d'optimisation serait bienvenue...
	/**
	 * Calcul de la dérivée première des polynomes de Lagrange
	 * \f[
	 *		\dot\phi_j(d_i)=\sum_{0\leq k\leq s,\;k\neq j}\left(\frac{1}{\alpha_j-\alpha_k}\times\prod_{0\leq l\leq s,\;l\neq k,\;l\neq j}\frac{d_i-\alpha_l}{\alpha_j-\alpha_l}\right)
	 * \f]
	 * Le résultat est un vecteur de taille dates.size() contenant pour chaque date
	 * un vecteur des valeurs des dérivées des s polynomes de Lagrange.
	 */
	void
	compute_polynomials_derivatives ()
	{
		int i; // index date
		int j; // index polynome
		int k; // index somme
		int l; // index produit
		//std::vector<std::vector<double>> ret;
		std::vector<double> tmp;
		std::vector<double> coeffs = chebyshev(m_degree);
		double somme;
		double prod;
		double t;

		m_polynomials_derivatives.clear();
		for (i=0; i<m_dates.size(); i++) {
			t = m_dates[i];
			tmp.clear();
			for (j=0; j<=m_degree; j++) {
				somme = 0.0;
				for (k=0; k<=m_degree; k++) {
					if (j!=k) {
						prod = 1.0;
						for (l=0; l<=m_degree; l++) {
							if ((l!=k) && (l!=j))
								prod *= (t-coeffs[l])/(coeffs[j]-coeffs[l]);
						}
						prod = prod/(coeffs[j]-coeffs[k]);
						somme += prod;
					}
				}
				tmp.push_back(somme);
			}
			m_polynomials_derivatives.push_back(tmp);
		}
	}
};

template <typename T_Q>
class LagrangeInterpolation : public LagrangePolynomials
{
public:
	std::vector<T_Q>
	pos_interp (int s, std::vector<double> dates, std::vector<T_Q> configs)
	{
		check_computed(s, dates);

		int i,j;
		std::vector<T_Q> ret;
		T_Q tmp;

		for (j=0; j<this->m_dates.size(); j++) {
			tmp = T_Q::Zero();
			for (i=0;i<=this->m_degree;i++) {
				tmp += m_polynomials[j][i]*configs[i];
			}
			ret.push_back(tmp);
		}

		return ret;
	}

	std::vector<T_Q>
	vel_interp (int s, std::vector<double> dates, std::vector<T_Q> configs, double h)
	{
		check_computed(s, dates);

		int i,j;
		std::vector<T_Q> ret;
		T_Q tmp;

		for (j=0; j<this->m_dates.size(); j++) {
			tmp = T_Q::Zero();
			for (i=0;i<=this->m_degree;i++) {
				tmp += m_polynomials_derivatives[j][i]*configs[i];
			}
			ret.push_back(tmp/h);
		}

		return ret;
	}
};

#endif
