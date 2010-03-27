#include <iostream>
#include <cmath>
#include <ctime>

#include "theoretic.hpp"

using namespace std;


int main() {
  BrownSim S;
  
  // int BrownSim::Sim(double T=t1, const int P2=25, const int Terms=100,
  // 		  const int Rb=4, // records begin with 2^Rb+1 points. 
  // 		  const int Re = P2-5, // records end with 2^Re+1 points.
  // 		  const int Nseg = 26,
  // 		  const int nRQ,
  // 		  double * RQuantile
  // 		  ) 

  time_t seconds = time (NULL);

  SimPara Para;
  Para.T = 1;
  Para.Terms = 100;
  Para.P2 = 23;
  Para.Nseg = 27;
  Para.Rb = 4;
  Para.Re = 20;
  Para.Rseed = seconds;
  //Para.IsQ = false;
  S.Sim(Para);
}

