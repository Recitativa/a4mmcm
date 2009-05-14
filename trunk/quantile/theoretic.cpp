#include "theoretic.hpp"
#include <cmath>
#include <iostream>


double Quatile::g(double x, double T) {
  double result;
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  result = g_inn(x,T, sigma, alpha, w);
  gsl_integration_workspace_free (w);
  return result;
}

inline double \
Quatile::g_inn(double x, double T, double sigma, double alpha,		\
	       gsl_integration_workspace *w) {
  
  double result, error;
  double Ialpha = 1.0;
  
  double params[] = {x, T, sigma, alpha};
  
  gsl_function F;
  
  std::cout << intg(1,params) << std::endl;
  F.function = &intg;
  F.params = params;
  
  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
			w, &result, &error); 
  return result;
}
// g_1(x-y;\alpha*t)*g_2(y;(1-\alpha)*t)
inline double Quatile::g1(double x,double t, double sigma) {
  if(x > 0)
    return 1/sigma * sqrt(2/PI*t)*exp(-x*x/2/(sigma*sigma)/t);
  else 
    return 0;
}
inline double Quatile::g2(double x, double t, double sigma) {
  return g1(-x,t, sigma);
}
inline double Quatile::intg(double y, void * Iparams) {
  double * params = (double *) Iparams;
  double x = params[0];
  double T = params[1];
  double sigma = params[2];
  double alpha = params[3];
  std::cout <<x <<T << sigma << alpha << std::endl;
  double f = g1(x-y,alpha*T, sigma)*g2(x, (1-alpha)*T, sigma);
  std::cout << f << std::endl;
  return f;
}


