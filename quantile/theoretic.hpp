


// Ref: Angelos Dassios, The Distribution of the quantile of a
// brownian motion with drift and the pricing of related
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

// path-dependent options, The annals of Applied Probability, 1995,
// Vol. 5, No. 2, 389-298

const double PI = 3.14159265358979323846;

class Brownian
{
public:
  double sigma;
  double mu;
}

class Quatile 
{
public:
  double alpha;
  Quatile(double Isigma = 1, double Ialpha=1) 
    : sigma(Isigma), alpha (Ialpha) {}
  double g(double x, double T);
  static double g_inn(double x, double T, double sigma, double alpha, \
		      gsl_integration_workspace *w);
private:
      // g_1(x-y;\alpha*t)*g_2(y;(1-\alpha)*t)
  static double g1(double x,double t, double sigma);
  static double g2(double x, double t, double sigma); 
  static double intg(double y, void * Iparams);
};

