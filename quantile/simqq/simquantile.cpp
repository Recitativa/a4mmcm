#include <iostream>
#include <cmath>
#include <ctime>
#include <sstream>

#include "theoretic.hpp"

using namespace std;


int main(int argc,      
          char *argv[],  
          char *envp[] )
{
  const int Narg = 8;

  if(argc != Narg) {
    cerr << argv[0] << " Terms Rb Rm Re sigma mu OutputIndex" << endl; 
    return 1;
  }
  stringstream ss (stringstream::in | stringstream::out);
  for(int i=1;i<Narg;i++) { ss << argv[i] << " ";}


  SimPara Para; // The class for Parameters 
  Para.T = 1; // Brownian motion from 0 to time T.
  ss >> Para.Terms; // simulate how many times.  
  // begin of the segements.
  ss >> Para.Rb;
  // middle of the segments.
  ss >> Para.Rm;
  // maximal segements, i.e. the end one. 
  ss >> Para.Re;
  // Seed for random number generator
  double sigma, mu;
  ss >> sigma >> mu;
  ss >> Para.Nseg; // index of the name of outputfile
  BrownSim S(sigma,mu);

  time_t seconds = time (NULL); // Get time as random seed. 
  Para.Rseed = seconds*Para.Nseg;

  S.Sim(Para);
}


