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
  BrownSim S;

  if(argc != 4) {
    cerr << argv[0] << " Terms Rb Re" << endl; 
    return 1;
  }
  stringstream ss (stringstream::in | stringstream::out);
  for(int i=1;i<9;i++) { ss << argv[i] << " ";}

  time_t seconds = time (NULL); // Get time as random seed. 

  SimPara Para; // The class for Parameters 
  Para.T = 1; // Brownian motion from 0 to time T.
  ss >> Para.Terms; // simulate how many times.  
  Para.P2 = 23; // Never use, for compatible.
  Para.Nseg = 27; // Never use, for compatible.
  // begin of the segements.
  ss >> Para.Rb;
  // maximal segements, i.e. the end one. 
  ss >> Para.Re;
  // Seed for random number generator
  Para.Rseed = seconds;
  S.Sim(Para);
}


