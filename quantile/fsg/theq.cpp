#include <cmath>
#include <iostream>
#include <gsl/gsl_integration.h>
     

using namespace std;

const double PI = 3.1415926535897932384626433832795;


double dP_dy(double t, double u, double y) {
  double s=y;
  return .5*t*((sqrt(2/PI/s)*exp(-u*u*s/2)		\
		-2*u*gsl_cdf_gaussian_Q(u*sqrt(s)))	\
	       *(2*u					\
		 + sqt(2/PI/(t-s))*exp(-u*u/2*(t-s))	\
		 -2*u*gsl_cdf_gaussian_Q(u*sqrt(s-t)))	\
	       );
}

double mu;
double sigma; //
double r; // interest rate
double S0; // initial value of stocks
double T; // time 
double K; // strike price

double f(double w, void *Iparams) { 
  return (S0-K)*exp(w)*dP_dy(T,mu/sigma,w);
}

double Quantile::mean() {

  gsl_integration_workspace * w

  = gsl_integration_workspace_alloc (3000);


  double result, error;

  double params[] = {};

  gsl_function F;

  F.function = &f;

  F.params = params;

  double alpha

  mu = r- sigma*sigma/2


  gsl_integration_qagi (&F, 1e-7, 1e-5, 1000,

                        w, &result, &error);

  gsl_integration_workspace_free(w);



  return result;

}



