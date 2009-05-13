#include "theoretic.hpp"
#include <cmath>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

// Ref: Angelos Dassios, The Distribution of the quantile of a
// brownian motion with drift and the pricing of related
// path-dependent options, The annals of Applied Probability, 1995,
// Vol. 5, No. 2, 389-298

class Quatile 
{
public:
  double sigma;
  double alpha;
  double T;
  Quatile(double Isigma = 1, double Ialpha=1, double IT=1) 
    : sigma(Isigma), alpha (Ialpha), T (IT) {}
  double g(double x) {
    gsl_integration_workspace * w 
      = gsl_integration_workspace_alloc (1000);
    
    double result, error;
    double expected = -4.0;
    double Ialpha = 1.0;
    
    gsl_function F;
    F.function = &p;
    F.params = &Ialpha;
    
    gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
			  w, &result, &error); 
  }
  
private:
  // g_1(x-y;\alpha*t)*g_2(y;(1-\alpha)*t)
  double g1(double x) {
    if(x > 0)
      return 1/sigma * sqrt(2/PI*T)*exp(-x*x/2/(sigma*sigma)/t);
    else 
      return 0;
  }
  double g2(double x) {
    return g1(-x);
  }
}
