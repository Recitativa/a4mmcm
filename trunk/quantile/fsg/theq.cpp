#include <cmath>
#include <iostream>
#include <gsl/gsl_integration.h>
     

using namespace std;

const double PI = 3.1415926535897932384626433832795;
double mu;
double sigma; //
double r; // interest rate
double S0; // initial value of stocks
double T; // time 
double K; // strike price



double dP_dy(double t, double u, double y) {
  double s=y;
  return .5*t*((sqrt(2/PI/s)*exp(-u*u*s/2)		\
		-2*u*gsl_cdf_gaussian_Q(u*sqrt(s)))	\
	       *(2*u					\
		 + sqt(2/PI/(t-s))*exp(-u*u/2*(t-s))	\
		 -2*u*gsl_cdf_gaussian_Q(u*sqrt(s-t)))	\
	       );
}

double phi(double t,x, u) {
}


gsl_integration_workspace * w			\
= gsl_integration_workspace_alloc (3000);

gsl_function iF;


int main() {
  double S0=100; // spot price
  double X=95;  //strike price
  double r=0.05; //interest rate
  double B; //down barrier
  double alpha=0.8; 
  double sigma=0.2;  //volatility
  double T=0.25;  //time to maturity 
  
  
  
}

double Quantile::mean() {



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



