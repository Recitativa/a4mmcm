#include "theoretic.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;

// Quantile, compute theoretical distribution for qunatiles 
// Never use. 
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




// Number of points recorded for the Brownian motion is 2<<Re
int BrownSim::Sim(SimPara Para) {
  // Type of real number, it should be long double if Re>20 since double 
  // has 53 bits for base, and approsimate 53/3 = 17 decimal digits, 
  // When Re>20, it need more than 17 digits in computation  
  typedef long double Real; 
  
  // Passing Parameters to local variables.  
  const double T = Para.T;
  const int Terms = Para.Terms;
  const int Rb = Para.Rb;
  const int Re = Para.Re;
  const int P2 = Para.P2;
  const int Nseg = Para.Nseg;
  const unsigned long int Rseed = Para.Rseed;


  int i; 
  
  int n= 1<<Re; // total number of segements
  double hsigma = sqrt(T/n); // corresponding sigma for each step; 

 
  // FileName format out_Rb_Re_.bin, a binary file. 
  ostringstream SoutFilename;
  SoutFilename << "nout_" << Rb << "_" << Re << "_";
  SoutFilename << ".bin";
  
  string outFilename = SoutFilename.str();
  ifstream testf;
  ofstream fout;
  testf.open(outFilename.c_str());

  
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
    int tRe = Re-1;
    fout.write((char *)&tRe, sizeof(int));    
    fout.write((char *)&P2, sizeof(int));    
    fout.write((char *)&Nseg, sizeof(int));
    fout.flush();
    cerr << "Output file " << outFilename << ": Head has written" << endl;
  } else {
    testf.close();
    fout.open(outFilename.c_str(), ios::app| ios::binary);
  }

  //setup random number generator. ref: GNU Scientific Library
  const gsl_rng_type * rngT;
  gsl_rng * r;
  rngT = gsl_rng_ran3;
  r = gsl_rng_alloc (rngT);
  gsl_rng_set(r, Rseed);

  
  Real B;
  double Q; // temporary varible for Quantile
  
  // the array stor the Brownina path
  Real * Record = new Real[1<<Re]; 
  
  if(Record == NULL) {
    cerr << "no enough memory!" << endl;
    return 1;
  }

  // Simulation
  for(int l = 0 ; l< Terms; l++) {
    B = 0; 
 
    // generate the path
    for(i=0; i< (1<< Re); i++) {
      B += (Real)gsl_ran_gaussian(r, hsigma);
      Record[i] = B;
    } 
    
    cerr << "coumputing Q" << endl;
    int nQ;
    int k;

    // compute Qunatiles from 1<<(Rb-1)/1<<Rb  to 1 step 1/1<<Rb
    // 
    // Consider Record as path with 1<<g points
    // e.g, g = Re-1, consider points labeled by "*"
    //      g = Re-2, consider points labeled by "o"
    // index means the index for the array. 
    // 
    // index 0    1    2    3                                     1<<Re-1
    //  0----|----|----|----|----|----|----|----|----|----|----|----|
    //            *         *         *         *         *         *     
    //                      o                   o                   o
    // as Np = 1<<g +1 points path
    for(int g=Rb; g<= Re; g++) {
      int Np = 1<<g;
      // following the example, nStep=2, when g = Re-1; 
      // nStep=4, when g= Re-2;
      size_t nStep = 1<<(Re-g);
      // StepIter is a Random Access Iterator which help STL 
      // consider Record as array with step nStep 
      StepIter<Real> Sp(Record+(nStep-1), nStep);
      // compute Quantiles for alpha: 1/2= 1<<(Rb-1)/1<<Rb, 
      //                                   1<<(Rb-1)/1<<Rb + 1/1<<Rb,
      //                                   .......
      //                              1  =  1<<Rb / 1<<Rb.
      // nQ is the index of Quantile when consider Record as a 1<<g path
      //       0    1     ...1<<(g-1)-1     1<<(g-1)-1+1<<(g-Rb)       1<<g-1
      //  0----|----|----|...----|----|----|----|----|----|----|----|----|
      //                         *                   *                   *
      //                         nQ                 nQ                  nQ
     for(k=0, nQ = 1<< (g-1);
	  k<= 1<<(Rb-1); k++, nQ += (1<< (g-Rb))) {
	nth_element (Sp, Sp+(nQ-1), Sp+Np);
	Q = (double)(Sp[nQ-1]);
	if(nQ == Np) Q = max(0.,Q);
	fout.write((char *)&Q, sizeof(double));	
      }
    }

    fout.flush();
    cerr << "Adding " << l << "th records" << endl;
  }
  fout.close();
  delete Record;
  gsl_rng_free (r);
  return 0;
}

