


// Ref: Angelos Dassios, The Distribution of the quantile of a
// brownian motion with drift and the pricing of related
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

#include <iostream>
using namespace std;

// path-dependent options, The annals of Applied Probability, 1995,
// Vol. 5, No. 2, 389-298

const double PI = 3.14159265358979323846;

class Brownian
{
public:
  double sigma;
  double mu;
  Brownian(double Isigma, double Imu): sigma(Isigma), mu(Imu) {};
};

class Quantile : public Brownian
{
public:
  double alpha;
  Quantile(double Isigma = 1, double Ialpha=1) 
    : Brownian(Isigma,0), alpha(Ialpha) {}
  double dX(double x);
  double mean();
private:
  static double dX_i(double x, void * Iparams);
  static double XdX_i(double x, void * Iparams) 
  { return x*dX_i(x, Iparams); }
};

