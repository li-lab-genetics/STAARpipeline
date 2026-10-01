// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

static double sparse_col_dot(const arma::sp_mat& A, const arma::sp_mat& B, arma::uword col)
{
	arma::sp_mat::const_col_iterator a_it = A.begin_col(col);
	arma::sp_mat::const_col_iterator a_end = A.end_col(col);
	arma::sp_mat::const_col_iterator b_it = B.begin_col(col);
	arma::sp_mat::const_col_iterator b_end = B.end_col(col);

	double out = 0;

	while ((a_it != a_end) && (b_it != b_end))
	{
		if (a_it.row() == b_it.row())
		{
			out += (*a_it) * (*b_it);
			++a_it;
			++b_it;
		}
		else if (a_it.row() < b_it.row())
		{
			++a_it;
		}
		else
		{
			++b_it;
		}
	}

	return out;
}

// [[Rcpp::export]]
List Individual_Score_Test_sp(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::vec residuals)
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

	int q = Sigma_iX.n_cols;

	arma::sp_mat Sigma_i_G;
	Sigma_i_G = Sigma_i*G;

	arma::mat tSigma_iX_G;
	tSigma_iX_G.zeros(q,p);
	tSigma_iX_G = trans(Sigma_iX)*G;

	for(i = 0; i < p; i++)
	{
		double Cov_ii = sparse_col_dot(G, Sigma_i_G, i);

		if (q > 0)
		{
			arma::vec tSigma_iX_G_i = tSigma_iX_G.col(i);
			Cov_ii = Cov_ii - arma::as_scalar(trans(tSigma_iX_G_i)*cov*tSigma_iX_G_i);
		}

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

