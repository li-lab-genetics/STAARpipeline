// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

static double sparse_dense_col_dot(const arma::sp_mat& A, const arma::mat& B, arma::uword col_a, arma::uword col_b)
{
	arma::sp_mat::const_col_iterator a_it = A.begin_col(col_a);
	arma::sp_mat::const_col_iterator a_end = A.end_col(col_a);

	double out = 0;

	while (a_it != a_end)
	{
		out += (*a_it) * B(a_it.row(), col_b);
		++a_it;
	}

	return out;
}

// [[Rcpp::export]]
List Individual_Score_Test_sp_denseGRM(arma::sp_mat G, const arma::mat& P, arma::vec residuals)
{
	int i;

	// number of markers
	int p = G.n_cols;

	// Uscore
	arma::rowvec Uscore = trans(residuals)*G;
	// log(p-value)
	arma::vec pvalue_log;
	pvalue_log.zeros(p);

	// SE of Uscore
	arma::vec Uscore_se;
	Uscore_se.zeros(p);
	// Effect size estimation
	arma::vec Est;
	Est.zeros(p);
	// SE of Effect size estimation
	arma::vec Est_se;
	Est_se.zeros(p);

	double test_stat = 0;

	arma::mat P_G;
	P_G = P*G;

	for(i = 0; i < p; i++)
	{
		double Cov_ii = sparse_dense_col_dot(G, P_G, i, i);

		if (Cov_ii == 0)
		{
			pvalue_log(i) = 0;
			Uscore_se(i) = 0;
			Est(i) = 0;
			Est_se(i) = 0;
		}
		else
		{
			test_stat = pow(Uscore(i),2)/Cov_ii;
			pvalue_log(i) = -R::pchisq(test_stat,1,false,true);

			Uscore_se(i) = sqrt(Cov_ii);
			Est(i) = Uscore(i)/Cov_ii;
			Est_se(i) = 1/Uscore_se(i);
		}

	}

	return List::create(Named("Score") = trans(Uscore), Named("Score_se") = Uscore_se, Named("pvalue_log") = pvalue_log, Named("Est") = Est, Named("Est_se") = Est_se);
}

