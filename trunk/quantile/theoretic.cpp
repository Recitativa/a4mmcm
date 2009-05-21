#include "theoretic.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;

double Quantile::dX(double x) {
  double params[] =  {sigma, alpha};
  return dX_i(x, params);
}

double Quantile::dX_i(double x, void * Iparams) {
  double * params = (double *)Iparams;
  //double sigma = params[0];
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
int BrownSim::Sim(double T=1, const int P2=20, const int Terms=100) {
  typedef double Real;
  const int Rb = 4; // records begin with 2^Rb+1 points.  
  const int Re = 10; // records end with 2^Re+1 points.


  QRCounter<Real, int> C1(-10,10, 1<<11);
  int i;
  int n= 1<<P2; // total number of segements
  double hsigma = T*sigma/n; // corresponding sigma for each step; 

  double RQuantile[] = {.5, .6, .7, .8, .9, 1.0};
  const int nRQ = sizeof(RQuantile) / sizeof(double); // #of Quantiles
  
  // FileName format out_P2_Rb_Re_quantiles
  ostringstream SoutFilename;
  SoutFilename << "out_" << P2 << "_ " << Rb << "_" << Re;
  for(i= 0 ; i< nRQ ; i++)
    SoutFilename << "_" << (int)(RQuantile[i]*100);
  SoutFilename <<".bin";
  string outFilename = SoutFilename.str();
  ifstream testf;
  ofstream fout;
  testf.open(outFilename.c_str());
  // output number of recorded quantile, and eqch quantile at first
  // write.
  // Format 
  // Head: int: P2 Rb Re nRQ double: RQuantiles
  // Record: Quantiles for total, 
  //         Quantiles for 1<< Rb number of elements
  //         Quantiles for 1<< (Rb+1) number of elements
  //         ...........................................
  //         Quantiles for 1<< Re number of elements 

  if(!testf.is_open()) {
    testf.close();
    fout.open(outFilename.c_str(), ios::app| ios::binary);
    fout.write((char *)&P2, sizeof(int));
    fout.write((char *)&Rb, sizeof(int));
    fout.write((char *)&Re, sizeof(int));
    fout.write((char *)&nRQ, sizeof(int));
    fout.write((char *)RQuantile, sizeof(RQuantile));
    fout.flush();
    cerr << "Output file " << outFilename << ": Head has written" << endl;
  } else {
    testf.close();
    fout.open(outFilename.c_str(), ios::app| ios::binary);
  }

  //setup random number generator
  const gsl_rng_type * rngT;
  gsl_rng * r;
  rngT = gsl_rng_taus;
  r = gsl_rng_alloc (rngT);
  gsl_rng_set(r, 123);

  
  Real B;
  double Q;
  Real * Record = new Real[(1<<Re)+1];
  if(Record == NULL) {
    cerr << "no enough memory!" << endl;
    return 1;
  }
  for(int l = 0 ; l< Terms; l++) {
    B = 0; 
    C1.init();
    Record[0]= B;
    C1.Add(B);
    for(i=0; i< (1<< Re); i++) {
      for(int j=0; j< 1<< (P2-Re); j++) {
	B += (Real)gsl_ran_gaussian(r, hsigma);
	C1.Add(B);
      }
      Record[i+1] = B;
    } 
    for(int k=0; k< nRQ; k++) {
      Q = (double)C1.QuantileC(RQuantile[k]);
      fout.write((char *)&Q, sizeof(double));
    }
    cerr << "coumputed Q" << endl;
    // as Np = 1<<g +1 points path
    for(int g=Rb; g< Re+1; g++) {
      const int Np = (1<<g)+1;
      StepIter<Real> Sp(Record,1<<(Re-g));
      int A;
      for(int k=0; k< nRQ; k++) {
	A = max(0,min(Np-1,(int)floor(Np*RQuantile[k])));
	nth_element (Sp, Sp+A, Sp+Np);
	Q = (double)(Sp[A]);
	fout.write((char *)&Q, sizeof(double));	
	//cerr << "G:" << k << endl;
      }
      //cerr << "Seg:" << g << endl;
    }
    fout.flush();
    cerr << "Adding " << l << "th records" << endl;
  }
  fout.close();
  delete Record;
  gsl_rng_free (r);
  return 0;
}

