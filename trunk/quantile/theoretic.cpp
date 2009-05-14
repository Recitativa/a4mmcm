#include "theoretic.hpp"
#include <cmath>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <fstream>

using namespace std;

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

void BrownSim::Sim(double T=1, int n=1000) {
  QRCounter<double, int> C1(-1,1, 1<<20);
  int i;
  double hsigma = T*sigma/n; 

  ofstream fout("out.bin", ios::app| ios::binary);
  //setup random number generator
  const gsl_rng_type * rngT;
  gsl_rng * r;
  rngT = gsl_rng_taus;
  r = gsl_rng_alloc (rngT);
  gsl_rng_set(r, 123);

  
  double B;
  for(int l = 0 ; l< 100; l++) {
    B = 0; i=0;
    C1.init();
    //C1.PrintDest();
    do {
      C1.Add(B);
      B += gsl_ran_gaussian(r, hsigma);
    } while(++i<n);
    double Q5 = C1.QuantileC(.5);
    fout.write((char *)&Q5, sizeof(double));
    fout.flush();
  }
  fout.close();
  gsl_rng_free (r);
}

