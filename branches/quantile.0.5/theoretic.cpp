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
int BrownSim::Sim(SimPara Para) {
  typedef long double Real;
  const double T = Para.T;
  const int Terms = Para.Terms;
  const int Rb = Para.Rb;
  const int Re = Para.Re;
  const int P2 = Para.P2;
  const int Nseg = Para.Nseg;
  const unsigned long int Rseed = Para.Rseed;

  int i;
  int n= 1<<P2; // total number of segements
  double hsigma = sqrt(T/n); // corresponding sigma for each step; 

 
  // FileName format out_P2_Rb_Re_quantiles
  ostringstream SoutFilename;
  SoutFilename << "sout_" << Rb << "_" << Re << "_";
  SoutFilename << P2 << "_" << Nseg << ".bin";
  
  string outFilename = SoutFilename.str();
  ifstream testf;
  ofstream fout;
  testf.open(outFilename.c_str());

  SoutFilename << ".txt";
  ofstream ftout;
  ftout.open(SoutFilename.str().c_str(), ios::out);

  // output number of recorded quantile, and eqch quantile at first
  // write.
  // Format 
  // Head: int: Rb Re P2 Nseg
  // Record: 
  //         Quantiles for 1<< Rb number of elements
  //         Quantiles for 1<< (Rb+1) number of elements
  //         ...........................................
  //         Quantiles for 1<< Re number of elements   
  //         Quantiles for 1<< P2 number of elements, i.e. the Dense one 

  if(!testf.is_open()) {
    testf.close();
    fout.open(outFilename.c_str(), ios::app| ios::binary);
    fout.write((char *)&Rb, sizeof(int));
    fout.write((char *)&Re, sizeof(int));    
    fout.write((char *)&P2, sizeof(int));    
    fout.write((char *)&Nseg, sizeof(int));
    fout.flush();
    cerr << "Output file " << outFilename << ": Head has written" << endl;
  } else {
    testf.close();
    fout.open(outFilename.c_str(), ios::app| ios::binary);
  }

  //setup random number generator
  const gsl_rng_type * rngT;
  gsl_rng * r;
  rngT = gsl_rng_ran3;
  r = gsl_rng_alloc (rngT);
  gsl_rng_set(r, Rseed);

  
  Real B;
  double Q;
  QRCounter<Real, int> C1(-10,10,1<<Nseg);
  
  Real * Record = new Real[1<<Re];
  if(Record == NULL) {
    cerr << "no enough memory!" << endl;
    return 1;
  }

  for(int l = 0 ; l< Terms; l++) {
    B = 0; 
    C1.init();
    for(i=0; i< (1<< Re); i++) {
      for(int j=0;j< (1 <<(P2-Re)); j++) {
	B += (Real)gsl_ran_gaussian(r, hsigma);
	C1.Add(B);
      }
      Record[i] = B;
    } 
    
    cerr << "coumputing Q" << endl;
    int nQ;
    int k;
    // Record from 1<<(Rb-1)/1<<Rb  to 1 step 1/1<<Rb
    
    // as Np = 1<<g +1 points path
    for(int g=Rb; g<= Re; g++) {
      int Np = 1<<g;
      size_t nStep = 1<<(Re-g);
      StepIter<Real> Sp(Record+(nStep-1), nStep);
      for(k=0, nQ = 1<< (g-1);
	  k<= 1<<(Rb-1); k++, nQ += (1<< (g-Rb))) {
	nth_element (Sp, Sp+(nQ-1), Sp+Np);
	Q = (double)(Sp[nQ-1]);
	//if(nQ == Np) Q = max(0.,Q);
	fout.write((char *)&Q, sizeof(double));	
	ftout << Q  << "\t";
	//cerr << "G:" << k << endl;
	//cerr << Q << " ";
      }
      ftout << endl;
      //cerr << "Seg:" << g << endl;
    }
    cerr << "computing Dense Quantiles.\t";
    for(k=0, nQ = 1<< (P2-1);
        k<= 1<<(Rb-1); k++, nQ += (1 <<(P2-Rb))) {
      try {
        Q = (double)C1.nQuantileC(nQ);
	//if(nQ == (1<<P2)) {
	//  Q = max(0., Q);
	//  cerr << "maximal in nQ:" << endl;
	//}
        fout.write((char *)&Q, sizeof(double));
        //cerr << Q << " ";
      } catch(OutofRangeException o) {
        cerr << "out of range when qantile: " << ((double)nQ)/(1<< P2);
      }
    }
    cerr << "coumputed Q" << endl;

    fout.flush();
    cerr << "Adding " << l << "th records" << endl;
  }
  fout.close();
  delete Record;
  gsl_rng_free (r);
  return 0;
}
