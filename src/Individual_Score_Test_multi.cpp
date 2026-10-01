// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Individual_Score_Test_multi(arma::mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::vec residuals, int n_pheno=1)
{
	int i,k,l;
	int p = G.n_cols;

	// number of markers
	int pp = p/n_pheno;

	// Uscore
	arma::rowvec Uscore = trans(residuals)*G;
	// log(p-value)
	arma::vec pvalue_log;
	pvalue_log.zeros(pp);

	arma::mat Uscore_cov;
	Uscore_cov.zeros(n_pheno,n_pheno);

	arma::uvec id_single;
	id_single.zeros(n_pheno);

	double test_stat = 0;

	int q = Sigma_iX.n_cols;

	arma::mat tSigma_iX_G;
	tSigma_iX_G.zeros(q,p);

	arma::mat tG_Sigma_i;
	tG_Sigma_i = trans(G)*Sigma_i;

	tSigma_iX_G = trans(Sigma_iX)*G;

	arma::mat quad;
	quad.zeros(1,1);

	for(i = 0; i < pp; i++)
	{
		for(k = 0; k < n_pheno; k++)
		{
			id_single(k) = k*pp+i;
		}

		for(k = 0; k < n_pheno; k++)
		{
			for(l = 0; l < n_pheno; l++)
			{
				Uscore_cov(k,l) = arma::as_scalar(tG_Sigma_i.row(id_single(k))*G.col(id_single(l)));
				if(q > 0)
				{
					Uscore_cov(k,l) -= arma::as_scalar(trans(tSigma_iX_G.col(id_single(k)))*cov*tSigma_iX_G.col(id_single(l)));
				}
			}
		}

		if (arma::det(Uscore_cov) == 0)
		{
			pvalue_log(i) = 0;
		}
		else
		{
			quad = trans(Uscore(id_single))*inv(Uscore_cov)*Uscore(id_single);
			test_stat = quad(0,0);
			pvalue_log(i) = -R::pchisq(test_stat,n_pheno,false,true);
		}

	}

	return List::create(Named("Score") = trans(Uscore), Named("pvalue_log") = pvalue_log);
}

