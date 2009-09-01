#include <iostream>
#include <fstream>
#include <cmath>     //standard mathematical library
#include <algorithm>
#include <vector>
#include <cassert>
using namespace std;

double factor=0.5;

double forwardST( double S, //spot price
		  double r, //interest rate
		  double B, //down barrier
		  double alpha, //
		  double sigma,  //volatility	
		  double T,  //time to maturity 
		  int steps) {
  double dt = T/steps;
  double R = exp((-r)*dt); 	//discount rate at each step
  double dw = sigma*dt; // up movement
  double du = (r)*dt; // up per unit time
  double D = alpha*(steps);  //barrier
  double pu = .5;  //probability of upward
  double pd = .5;   //downward
  double p0 = 0;
  int m = steps;
  double bar = log(B/S); //barrier
  const int N = (m+1)*(2*m+1);
  double data[2*N];

#define g(k,j) ((dw*(j)+du*i<=bar)?((k)+1):(k))
#define IND(k,j)  ((k)*(2*m+1)+(j)+m)
	
  double * even, * odd, *tmp;
  even  = data;
  odd = data+N;

  for(int k=0; k<m+1;k++)
    for(int j=-m; j<m+1;j++) {
      if (k<D) even[IND(k,j)]=1;
      else even[IND(k,j)]=0;
    }
	
  for(int i=m-1;i >=0; i--) {
    tmp = even; even = odd; odd = tmp;
    for(int k=0; k<i+1;k++)
	for(int j=-i; j<i+1;j++) {
	  even[IND(k,j)] = (  pu*odd[IND(g(k,j+1),j+1)]			\
			      + p0*odd[IND(g(k,j),j)]			\
			      + pd*odd[IND(g(k,j-1),j-1)])*R;
	}
  }  
  return even[IND(0,0)];
}

int main() {
  double S=100; //spot price
  double X=95;  //strike price
  double r=0.05; //interest rate
  double B; //down barrier
  double alpha=0.8; 
  double sigma=0.2;  //volatility
  double T=0.25;  //time to maturity 
  int steps=10;
  double value;

  ofstream fp;
  fp.open("fsg.dat");


  for( steps=10; steps<=100; steps+=10)
  {
    double valp, valn;
    valp = 0;
    valn = exp(-r*T);
    value = 0;
    int n = 100;
    double dd = ((r-sigma*sigma/2)*T+sigma*T)/n;
    for(int k=log(X/S)/dd; k<n; k++) {
      B= S*exp(k*dd);
      valp = valn;
      valn=forwardST(S, r, B, alpha, sigma, T, steps); 
      //cerr<<k <<" "<< B<< " "<< valp << " " << valn << endl;
      value=value+max(B-X,0.0)*(valp-valn);
    }
    fp << steps << "  " << value << endl;
  }		
  fp.close();
  return 0;
}
