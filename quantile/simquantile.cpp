#include <iostream>
#include <cmath>
#include <ctime>

#include "theoretic.hpp"

using namespace std;


int main() {
  BrownSim S;
  

  time_t seconds = time (NULL); // Get time as random seed. 

  SimPara Para; // The class for Parameters 
  Para.T = 1; // Brownian motion from 0 to time T.
  Para.Terms = 200; // simulate how many times.  
  Para.P2 = 23; // Never use, for compatible.
  Para.Nseg = 27; // Never use, for compatible.
  // begin of the segements.
  Para.Rb = 6;
  // maximal segements, i.e. the end one. 
  Para.Re = 20;
  // Seed for random number generator
  Para.Rseed = seconds;
  S.Sim(Para);
}


