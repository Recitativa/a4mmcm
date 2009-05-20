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

// P2 number of points in brwonian motion is 2<<P2
void BrownSim::Sim(double T=1, const int P2=20) {
  typedef double Real;
  const int Rb = 4; // records begin with 2^Rb+1 points.  
  const int Re = 21; // records end with 2^Re+1 points.

  Real * Record = new real[1<<Re+1];

  QRCounter<Real, int> C1(-10,10, 1<<20);
  int i;
  int n= 2<<P2;
  double hsigma = T*sigma/n; 

  double RQuantile[] = {.5, .6, .7, .8, .9, 1.0};
  const int nRQ = sizeof(RQuantile) / sizeof(double);
  
  ofstream fout("out.bin", ios::app| ios::binary);
  //write number of recorded quantile, and eqch quantile
  fout.write((char *)&nRQ, sizeof(int));
  fout.flush((char *)RQuantile, sizeof(RQuantile));

  //setup random number generator
  const gsl_rng_type * rngT;
  gsl_rng * r;
  rngT = gsl_rng_taus;
  r = gsl_rng_alloc (rngT);
  gsl_rng_set(r, 123);

  
  Real B;
  for(int l = 0 ; l<100; l++) {
    B = 0; i=0;
    C1.init();
    Record[0]= B;
    C1.Add(B);
    for(i=0; i< 1<< Re; i++) {
      for(j=0; j< 1<< (P2-Re); j++) {
	B += (Real)gsl_ran_gaussian(r, hsigma);
	C1.Add(B);
      }
      Record[i+1] = B;
    } 
    Real Q5 = C1.QuantileC(.5);
    fout.write((char *)&Q5, sizeof(double));
    fout.flush();
  }
  fout.close();
  gsl_rng_free (r);
}

