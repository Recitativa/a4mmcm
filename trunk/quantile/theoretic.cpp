#include "theoretic.hpp"
#include <cmath>
#include <iostream>




double Quantile::dX(double x) {
  double params[] =  {sigma, alpha};
  return dX_i(x, params);
}

double Quantile::dX_i(double x, void * Iparams) {
  double * params = (double *)Iparams;
  double sigma = params[0];
  double alpha = params[1];
  double beta = sqrt((1. - alpha) / alpha);
  static const double c1 = sqrt(2 / PI) * 2;

  if (x >= 0)
    return c1 * exp(-x*x / 2) * gsl_cdf_ugaussian_Q(beta * x);
  else
    return c1 * exp(-x*x / 2) * gsl_cdf_ugaussian_Q(-x / beta);
}

double Quantile::mean() {
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  double result, error;
  double params[] = {sigma, alpha};

  gsl_function F;
  F.function = &XdX_i;
  F.params = params;

  gsl_integration_qagi (&F, 1e-7, 1e-5, 1000,
                        w, &result, &error);
  gsl_integration_workspace_free(w);

  return result;
}

template<class Num, class C>
QRCounter<Num, C>::QRCounter(Num begin, Num end, size_t In): Rbegin(begin), Rend(end), n(In), total(0) {
  counter = new C[n](0);
}

template<class Num, class C>
QRCounter<Num, C>::~QRCounter() {
  delete counter;
}

template<class Num, class C>
int QRCounter<Num, C>::Add(Num x) {
  return 0;
}

template<class Num, class C>
void QRCounter<Num, C>::PrintDest() {
  for(size_t i =0; i< n; i++) 
    std::cout << counter[i] << "\t";
}

