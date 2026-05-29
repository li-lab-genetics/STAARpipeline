// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// declare K_Binary_SPA_wt
double K_Binary_SPA_wt(double x, arma::vec muhat, arma::vec G);
// declare K_Binary_SPA_alt_wt (alternative way if not converge)
double K_Binary_SPA_alt_wt(double x, arma::vec muhat, arma::vec G);
// declare K1_Binary_SPA_wt (first derivative)
double K1_Binary_SPA_wt(double x, arma::vec muhat, arma::vec G, double q);
// declare K1_Binary_SPA_alt_wt (first derivative, alternative way if not converge)
double K1_Binary_SPA_alt_wt(double x, arma::vec muhat, arma::vec G, double q);
// declare K2_Binary_SPA_wt (second derivative)
double K2_Binary_SPA_wt(double x, arma::vec muhat, arma::vec G);
// declare K2_Binary_SPA_alt_wt (second derivative, alternative way if not converge)
double K2_Binary_SPA_alt_wt(double x, arma::vec muhat, arma::vec G);
// declare NR
double NR_Binary_SPA_wt(arma::vec muhat, arma::vec G, double q, double init, double tol, int max_iter);
// declare golden section search method for selecting initial value of bisection algorothm
arma::vec goldenSectionSearchForSignChange_wt(double a, double b, arma::vec muhat, arma::vec G, double q, double tol, int max_iter);
// declare Bisection method
double Bisection_Binary_SPA_wt(arma::vec muhat, arma::vec G, double q, double xmin, double xmax, double tol);
// declare Saddle_Binary_SPA_wt using NR
double Saddle_Binary_SPA_wt(double q, arma::vec muhat, arma::vec G, double tol, int max_iter, bool lower);
// declare Saddle_Binary_SPA_wt using Bisection
double Saddle_Binary_SPA_Bisection_wt(double q, arma::vec muhat, arma::vec G, double tol, int max_iter, double xmin, double xmax, bool lower);
// check na
bool check_is_na_wt(double x); 
// check same sign
bool haveSameSign_wt(double a, double b);

// [[Rcpp::export]]
arma::vec Individual_Score_Test_SPA_wt(arma::mat G, arma::vec residuals, arma::vec muhat, double tol, int max_iter)
{
	int i,k;
	double sum0 = 0.0;
	
	bool lower1 = false;
	bool lower2 = true;

	int un = G.n_cols;

	// p-values
	arma::vec res;
	res.ones(un);
	
	// p-values components
	double respart1 = 0.0;
	double respart2 = 0.0;
	
	// Score statistics
	arma::rowvec x = trans(residuals)*G;
	int n = x.size();
	
	// init value of bisection
	double xmin = -100.0;
	double xmax = 100.0;

	for(i = 0; i < un; i++)
	{
		sum0 = x(i);
		
		// calculate p-value
		respart1 = 1.0;
		respart2 = 1.0;
		respart1 = Saddle_Binary_SPA_wt(fabs(sum0), muhat, G.col(i), tol, max_iter, lower1);
		
		if((check_is_na_wt(respart1))||(respart1 == 1.0))
		{
			respart1 = Saddle_Binary_SPA_Bisection_wt(fabs(sum0), muhat, G.col(i), tol, max_iter, xmin, xmax, lower1);
		}
		
		if((check_is_na_wt(respart1))||(respart1 == 1.0))
		{
			res(i) = 1;
		}else
		{
			respart2 = Saddle_Binary_SPA_wt(-fabs(sum0), muhat, G.col(i), tol, max_iter, lower2);
			
			if((check_is_na_wt(respart2))||(respart2 == 1.0))
			{
				respart2 = Saddle_Binary_SPA_Bisection_wt(-fabs(sum0), muhat, G.col(i), tol, max_iter, xmin, xmax, lower2);
			}
			
			if((check_is_na_wt(respart2))||(respart2 == 1.0))
			{
				res(i) = 1;
			}else
			{
				res(i) = respart1 + respart2;
			}
		}
		
		// if saddle point approximation fails, set the p-value as 1.0
		if(res(i) > 1)
		{
			res(i) = 1;		
		}
	}

	return res;
}


double Saddle_Binary_SPA_wt(double q, arma::vec muhat, arma::vec G, double tol, int max_iter, bool lower)
{
	bool logp = false;
	
	// init for NR
	double init = 0.0;

	double ki = 0.0;
	double res = 0.0;

	// Saddle point approximation
    double xhat = NR_Binary_SPA_wt(muhat,G,q,init,tol,max_iter);
    double w = sqrt(2*(xhat*q - K_Binary_SPA_wt(xhat,muhat,G)));
	
	if((R_finite(w)==0)||(check_is_na_wt(w)))
	{
		w = sqrt(2*(xhat*q - K_Binary_SPA_alt_wt(xhat,muhat,G)));
	}
	
    if (xhat < 0){
        w = -w;
    }

	ki = xhat*sqrt(K2_Binary_SPA_wt(xhat,muhat,G));
	if((R_finite(ki)==0)||(check_is_na_wt(ki)))
	{
		ki = xhat*sqrt(K2_Binary_SPA_alt_wt(xhat,muhat,G));
	}
	
	if(fabs(xhat)<1e-04)
	{
		res = 1.0;
	}else
	{
		res =  R::pnorm(w+log(ki/w)/w,0.0,1.0,lower,logp);
	}
    
	return res;
}


double NR_Binary_SPA_wt(arma::vec muhat, arma::vec G, double q, double init, double tol, int max_iter)
{
	double xi = 0.0;
	double xi_update = 0.0;
	
	// first step
	xi = init;
	xi_update = init;
	
	if(fabs(K1_Binary_SPA_wt(xi, muhat, G, q)) > tol)
	{
		xi_update = xi - K1_Binary_SPA_wt(xi, muhat, G, q)/K2_Binary_SPA_wt(xi, muhat, G);
	}
	
	// iteration number
	int no_iter = 0;
	double numerator = 0.0;
	double denominator = 1.0;
	
	while (!(check_is_na_wt(xi_update))&&(R_finite(xi_update)==1)&&(fabs(xi_update - xi) > tol)&&(fabs(K1_Binary_SPA_wt(xi_update, muhat, G, q)) > tol)&&(no_iter < max_iter))
    {
		no_iter = no_iter + 1;
		xi = xi_update;	
		
		// calculate numerator
		numerator = K1_Binary_SPA_wt(xi, muhat, G, q);
		if((R_finite(numerator)==0)||(check_is_na_wt(numerator)))
		{
			numerator = K1_Binary_SPA_alt_wt(xi, muhat, G, q);
		}
		
		// calculate denominator
		denominator = K2_Binary_SPA_wt(xi, muhat, G);
		if((R_finite(denominator)==0)||(check_is_na_wt(denominator)))
		{
			denominator = K2_Binary_SPA_alt_wt(xi, muhat, G);
		}
		
		xi_update = xi - numerator/denominator;
    }

	// test finite
	if((R_finite(xi_update)==0)||(check_is_na_wt(xi_update)))
	{
		xi_update = xi;
	}
	
	
	// test convergence
	double w = 0.0; 
	double xhat = 0.0; 
	double ki = 0.0;

	if(no_iter == max_iter)
	{
		w = sqrt(2*(xhat*q - K_Binary_SPA_wt(xhat,muhat,G)));
		xhat = xi_update;
		
		if (xhat < 0){
			w = -w;
		}

		ki = xhat*sqrt(K2_Binary_SPA_wt(xhat,muhat,G));
		if((R_finite(ki)==0)||(check_is_na_wt(ki)))
		{
			ki = xhat*sqrt(K2_Binary_SPA_alt_wt(xhat,muhat,G));
		}
	
		if(fabs(w+log(ki/w)/w) < 38)
		{
			xi_update = 0.0;	
		}
	}
	
	return xi_update;
}

double Saddle_Binary_SPA_Bisection_wt(double q, arma::vec muhat, arma::vec G, double tol, int max_iter, double xmin, double xmax, bool lower)
{
	bool logp = false;
	
	double ki = 0.0;
	double res = 1.0;

	// Saddle point approximation
	// golden section search for selecting init value of bisection algorithm
	arma::vec xlimit;
	xlimit.ones(2);
	
	xlimit = goldenSectionSearchForSignChange_wt(xmin, xmax, muhat, G, q, tol, max_iter);
	
	if(xlimit(0) < xlimit(1))
	{
		xmin = xlimit(0);
		xmax = xlimit(1);
	}else
	{
		xmin = xlimit(1);
		xmax = xlimit(0);
	}
	
    double xhat = Bisection_Binary_SPA_wt(muhat, G, q, xmin, xmax, tol);
    double w = sqrt(2*(xhat*q - K_Binary_SPA_wt(xhat,muhat,G)));
	
	if((R_finite(w)==0)||(check_is_na_wt(w)))
	{
		w = sqrt(2*(xhat*q - K_Binary_SPA_alt_wt(xhat,muhat,G)));
	}
	
    if (xhat < 0){
        w = -w;
    }

	ki = xhat*sqrt(K2_Binary_SPA_wt(xhat,muhat,G));
	if((R_finite(ki)==0)||(check_is_na_wt(ki)))
	{
		ki = xhat*sqrt(K2_Binary_SPA_alt_wt(xhat,muhat,G));
	}
	
	if(fabs(xhat)<1e-04)
	{
		res = 1.0;
	}else
	{
		res =  R::pnorm(w+log(ki/w)/w,0.0,1.0,lower,logp);
	}
    
	return res;
}



double Bisection_Binary_SPA_wt(arma::vec muhat, arma::vec G, double q, double xmin, double xmax, double tol)
{

    // the range of x to search
    double xupper = xmax;
    double xlower = xmin;
	
	double K1left = 0.0;
	double K1right = 0.0;
	
    double x0 = 0.0;
    double K1x0 = 1.0;
	
	K1left = K1_Binary_SPA_wt(xlower, muhat, G, q);
	if((R_finite(K1left)==0)||(check_is_na_wt(K1left)))
	{
		K1left = K1_Binary_SPA_alt_wt(xlower, muhat, G, q);
	}
	
	K1right = K1_Binary_SPA_wt(xupper, muhat, G, q);
	if((R_finite(K1right)==0)||(check_is_na_wt(K1right)))
	{
		K1right = K1_Binary_SPA_alt_wt(xupper, muhat, G, q);
	}
	
	if((!((R_finite(K1right)==0)||(check_is_na_wt(K1right))))&&(!((R_finite(K1left)==0)||(check_is_na_wt(K1left)))))
	{	
		if(!haveSameSign_wt(K1left, K1right))
		{
			while ((fabs(xupper-xlower) > tol)&&(fabs(K1x0) > tol))
			{
				x0 = (xupper + xlower)/2.0;
				K1x0 = K1_Binary_SPA_wt(x0, muhat, G, q);
				
				if((R_finite(K1x0)==0)||(check_is_na_wt(K1x0)))
				{
					K1x0 = K1_Binary_SPA_alt_wt(x0, muhat, G, q);
				}
			
				if(K1x0 == 0)
				{
					break;
				}else
				{
					if(haveSameSign_wt(K1left, K1x0))
					{
						xlower = x0;
					}else
					{
						xupper = x0;
					}
				
					K1left = K1_Binary_SPA_wt(xlower, muhat, G, q);
					if((R_finite(K1left)==0)||(check_is_na_wt(K1left)))
					{
						K1left = K1_Binary_SPA_alt_wt(xlower, muhat, G, q);
					}
	
					K1right = K1_Binary_SPA_wt(xupper, muhat, G, q);
					if((R_finite(K1right)==0)||(check_is_na_wt(K1right)))
					{
						K1right = K1_Binary_SPA_alt_wt(xupper, muhat, G, q);
					}
				}
			}
		}		
	}
    return x0;
}

arma::vec goldenSectionSearchForSignChange_wt(double a, double b, arma::vec muhat, arma::vec G, double q, double tol, int max_iter)
{
	int iter = 0;

    // results init
	arma::vec x;
	x.ones(2);
	x(0) = a;
	x(1) = b;

	// K1 init
	double K1x1 = 0.0;
	double K1x0 = 0.0;
	
	double phi = (1 + sqrt(5)) / 2;
	
	K1x0 = K1_Binary_SPA_wt(x(0), muhat, G, q);
	if((R_finite(K1x0)==0)||(check_is_na_wt(K1x0)))
	{
		K1x0 = K1_Binary_SPA_alt_wt(x(0), muhat, G, q);
	}
	
	K1x1 = K1_Binary_SPA_wt(x(1), muhat, G, q);
	if((R_finite(K1x1)==0)||(check_is_na_wt(K1x1)))
	{
		K1x1 = K1_Binary_SPA_alt_wt(x(1), muhat, G, q);
	}
	
	while ((iter < max_iter) && (fabs(x(1) - x(0)) > tol) && haveSameSign_wt(K1x0,K1x1)) 
	{
		x(1) = b - (b - a) / phi;
        x(0) = a + (b - a) / phi;
		
		K1x0 = K1_Binary_SPA_wt(x(0), muhat, G, q);
		if((R_finite(K1x0)==0)||(check_is_na_wt(K1x0)))
		{
			K1x0 = K1_Binary_SPA_alt_wt(x(0), muhat, G, q);
		}
	
		K1x1 = K1_Binary_SPA_wt(x(1), muhat, G, q);
		if((R_finite(K1x1)==0)||(check_is_na_wt(K1x1)))
		{
			K1x1 = K1_Binary_SPA_alt_wt(x(1), muhat, G, q);
		}
	
        if (K1x1 < K1x0) 
		{
            b = x(0);
        }else 
		{
            a = x(1);
        }
        iter++;
    }
	
	return x;
	
}



double K_Binary_SPA_wt(double x, arma::vec muhat, arma::vec G)
{
    double res = 0.0;
    const int n = muhat.size();

    for(int i = 0; i < n; i++)
    {
		res = res - x * muhat(i) * G(i);
        res = res + log(1 - muhat(i) + muhat(i) * exp(x * G(i)));
    }

    return res;
}



double K_Binary_SPA_alt_wt(double x, arma::vec muhat, arma::vec G)
{
    double res = 0.0;
    const int n = muhat.size();

    for(int i = 0; i < n; i++)
    {
		// res = res - x * muhat(i) * G(i);
		
        res = res + log((1 - muhat(i))*exp(-x * G(i)) + muhat(i));
    }

    return res;
}


double K1_Binary_SPA_wt(double x, arma::vec muhat, arma::vec G, double q)
{
    double res = 0.0;
    const int n = muhat.size();

    for(int i = 0; i < n; i++)
    {
        res = res - muhat(i) * G(i);
		res = res + muhat(i) * G(i)/(muhat(i) + (1 - muhat(i))*exp(-x * G(i)));
    }
	
    res = res - q;

    return res;
}



double K1_Binary_SPA_alt_wt(double x, arma::vec muhat, arma::vec G, double q)
{
    double res = 0.0;
    const int n = muhat.size();

    for(int i = 0; i < n; i++)
    {
        res = res - muhat(i) * G(i);
		res = res + muhat(i) * G(i) * exp(x * G(i)) /(muhat(i)* exp(x * G(i)) + (1 - muhat(i)));
    }
	
    res = res - q;

    return res;
}


double K2_Binary_SPA_wt(double x, arma::vec muhat, arma::vec G)
{
    double res = 0.0;
	double temp1 = 0.0;
	double temp2 = 1.0;
	
    const int n = muhat.size();

    for(int i = 0; i < n; i++)
    {
		temp1 = muhat(i) * (1 - muhat(i)) * pow(G(i),2.0)*exp(-x * G(i));
		temp2 = muhat(i) + (1 - muhat(i)) * exp(-x * G(i));
		
        res = res + temp1/pow(temp2,2.0);
    }
	
    return res;
}


double K2_Binary_SPA_alt_wt(double x, arma::vec muhat, arma::vec G)
{
    double res = 0.0;
	double temp1 = 0.0;
	double temp2 = 1.0;
	
    const int n = muhat.size();

    for(int i = 0; i < n; i++)
    {
		temp1 = muhat(i) * (1 - muhat(i)) * pow(G(i),2.0);
		temp2 = muhat(i) * exp(x * G(i)) + (1 - muhat(i));
		
        res = res + temp1/pow(temp2,2.0);
    }
	
    return res;
}

bool check_is_na_wt(double x) {
  return Rcpp::traits::is_na<REALSXP>(x);
}

bool haveSameSign_wt(double a, double b) {
    return (a >= 0) == (b >= 0);
}