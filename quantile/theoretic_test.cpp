#define BOOST_TEST_DYN_LINK
#define BOOST_AUTO_TEST_MAIN
#include <iostream>
#include <boost/test/auto_unit_test.hpp>
#include <cmath>
#include <sstream>
#include <fstream>

#include "theoretic.hpp"
#include <gsl/gsl_integration.h>

using namespace std;

double f(double x, void * params) {
  double alpha = *(double *) params;
  double ff = log(alpha * x) / sqrt(x);
  return ff;
}

double inte_test() {
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  double result, error;
  double expected = -4.0;
  double alpha = 1.0;

  gsl_function F;
  F.function = &f;
  F.params = &alpha;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                        w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}


BOOST_AUTO_TEST_CASE( test1 ) {
  // with or without some check
  BOOST_REQUIRE( fabs(inte_test() + 4.0) < 0.01);

  Quantile Q1(1.0, .3);

  double dX1_E = 0.0612817;
  double mean_E =  -0.230539;
  double dX1 = Q1.dX(1);
  double mean = Q1.mean();
  cout << "dX(1) = " << dX1_E << endl;
  cout << "mean " << mean_E << endl;

  BOOST_REQUIRE( fabs(mean_E - mean) < 0.001);
  BOOST_REQUIRE( fabs(dX1_E - dX1) < 0.001);

  QRCounter<double, long> C1(2,8,1<<20);
  double Data[12] = {2, 3, 4, 5, 1, 
		     6, 9, 11, 8, 7, 
		     12, 10};
  for(int i=0; i< 12; i++)
    C1.Add(Data[i]);
  //C1.PrintDest();
  cerr << "finish adding" << endl;
  try {
    cerr << "Countinous Quantile .5 : " << C1.QuantileC(.5) << endl;
    cerr << "Quantile .5 : " << C1.Quantile(.5) << endl;
    cerr << "Quantile .5 : " << C1.QuantileC(.9) << endl;
    BOOST_REQUIRE(false && "never go here");
  } 
  catch( OutofRangeException e ) {
    if (e.i ==0) 
      cerr << "Out of lower range" << endl;
    else
      cerr << "Out of upper range" << endl;
  }
  //QRCounter<double, long> C1(2,8,268435456);
  //BrownSim S1;
  //S1.Sim(1,100);

  StepIter<double> Sp1(Data,2), Sp2, Sp3(Data,3);  
  Sp2 = Sp1+1;
  cerr << "Sp1[0]: " << Sp1[0] << endl;
  BOOST_REQUIRE(Sp1[0]==2);
  cerr << "Sp2: " << Sp2[0]  << endl;
  BOOST_REQUIRE(Sp2[0]==4);
  cerr << "Sp2+1: " << *(Sp2+1)  << endl;  
  BOOST_REQUIRE(*(Sp2+1)==1);
  try {
    cerr << Sp3-Sp2;
    BOOST_REQUIRE(false && "never access here: Sp3-Sp2 for different 'steps'");
  } catch(DifferentStepsException) {
    cerr << "Different Steps Expected!" << endl;
  }
  Sp2 += 1;
  cerr << "Sp2:" << *Sp2 << endl;
  cerr << "Sp2-Sp1 " << Sp2-Sp1 << endl;
  Sp2 -=1;
  cerr << "Sp2:" << *Sp2 << endl;
  cerr << "Sp2-Sp1 " << Sp2-Sp1 << endl;
  Sp3 = Sp2-1;
  cerr << "Sp3:" << *Sp3 << endl;
  cerr << "Sp3-Sp2 " << Sp3-Sp2 << endl;
  cerr << "Sp2<Sp3 ? :" << (Sp2<Sp3) << endl;
  cerr << "Sp2>Sp3 ? :" << (Sp2>Sp3) << endl;
  cerr << "Sp3<Sp1 ? :" << (Sp3<Sp1) << endl;
  cerr << "Sp3<Sp1 ? :" << (Sp3<=Sp1) << endl;
  cerr << "*(Sp3++) :" << *(Sp3++) << endl;
  cerr << "*Sp3 :" << *Sp3 << endl;
  Sp1 = Data+1;
  cerr << "Sp1=Data+1 : " << *Sp1 << endl;
  
#define PRT_DATA do {				\
    for(int i=0; i<12; i++)			\
      cerr << Data[i] << " ";			\
    cerr << endl;				\
  } while(false)
  PRT_DATA;
  StepIter<double> Sp_4(Data, 4);
  sort(Sp_4, Sp_4+3);
  PRT_DATA;
  StepIter<double> Sp_3(Data, 3);
  sort(Sp_3, Sp_3+4);
  PRT_DATA;
  StepIter<double> Sp_2(Data, 2);
  sort(Sp_2, Sp_2+6);
  PRT_DATA;
  StepIter<double> Sp_1(Data, 1);
  sort(Sp_1, Sp_1+12);
  PRT_DATA;
  StepIter<double> Sp_r(Data+11, -1);
  sort(Sp_r, Sp_r+12);
  PRT_DATA;


  
  // open files
  int i;
  double RQuantile[] = {.5, .6, .7, .8, .9, 1.0};
  const int nRQ = sizeof(RQuantile) / sizeof(double);
  
  ostringstream SoutFilename;
  SoutFilename << "out";
  for(i= 0 ; i< nRQ ; i++)
    SoutFilename << "_" << (int)(RQuantile[i]*100);
  SoutFilename <<".bin";
  string outFilename = SoutFilename.str();
  ifstream testf(outFilename.c_str());
  ofstream fout(outFilename.c_str(), ios::app| ios::binary);

 
  // output number of recorded quantile, and eqch quantile at first
  // write.
  if(!testf.is_open()) {
    fout.write((char *)&nRQ, sizeof(int));
    fout.write((char *)RQuantile, sizeof(RQuantile));
    fout.flush();
  }
  testf.close();
  for(int l=0; l< 10000; l++) {
    double Q = 3.213123124325324324242424*l;
    fout.write((char *)&Q, sizeof(double));
    fout.flush();
  }
  fout.close();

#define PRT_DD do {				\
    for(int i=0; i<9; i++)			\
      cerr << DD[i] << " ";			\
    cerr << endl;				\
  } while(false)

  
  double DD[] = {8,3,1,9,2,5,4,6,7};
  StepIter<double> SpD2(DD,2), SpD1(DD,1);
  nth_element(SpD2, SpD2+2, SpD2+5); 
  PRT_DD;
  nth_element(SpD1, SpD1+4, SpD1+9); 
  PRT_DD;
 

}


