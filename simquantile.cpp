#include <iostream>
#include <cmath>

#include "theoretic.hpp"

using namespace std;


int main() {
  BrownSim S;
  double RQuantile[] = {8./16, 9./16, 10./16, 11./16,12./16,13./16,14./16,
			15./16, 16./16};
  const int nRQ = sizeof(RQuantile) / sizeof(double); // #of Quantiles

  // int BrownSim::Sim(double T=t1, const int P2=25, const int Terms=100,
  // 		  const int Rb=4, // records begin with 2^Rb+1 points. 
  // 		  const int Re = P2-5, // records end with 2^Re+1 points.
  // 		  const int Nseg = 26,
  // 		  const int nRQ,
  // 		  double * RQuantile
  // 		  ) 

  SimPara Para;
  Para.nRQ = nRQ;
  Para.RQuantile = RQuantile;
  Para.T = 1;
  Para.P2 = 22;
  Para.Terms = 100;
  Para.Rb = 4;
  Para.Re = Para.P2-5;
  Para.Nseg = 22;
  Para.Rseed = 123;
  S.Sim(Para);
}


