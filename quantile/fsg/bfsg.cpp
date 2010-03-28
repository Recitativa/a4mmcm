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
  double R = exp(-r*dt); 	//discount rate at each step
  double dw = sigma*sqrt(dt); // up movement
  double du = (r-sigma*sigma/2)*dt; // up per unit time
  double D = alpha*(steps+1);  //barrier
  double pu = .5;  //probability of upward
  double pd = .5;   //downward
  double p0 = 0;
  int m = steps;
  double bar = log(B/S); //barrier
  const int N = (m+1)*(2*m+1);

  double odd[m+1][2*m+1];
  double even[m+1][2*m+1];

#define g(k,j) ((dw*(j)+du*i<=bar)?((k)+1):(k))
	
  for(int k=0; k<m+1;k++) {
    for(int j=-m; j<m+1;j++) {
      if (k<D) even[k][j+m]=1;
      else even[k][j+m]=0;
    }
  }
	
  for(int i=m-1;i >=0; i--) {
    for(int k=0; k<i+1;k++)
      if((m-i)%2==1) {
	for(int j=-i; j<i+1;j++) {
	  odd[k][j+m]=(pu*even[g(k,j+1)][j+1+m]+p0*even[g(k,j)][j+m]+pd*even[g(k,j-1)][j-1+m])*R; }
      }
      else{ 
	for(int j=-i; j<i+1;j++) {
	  even[k][j+m]=(pu*odd[g(k,j+1)][j+1+m]+p0*odd[g(k,j)][j+m]+pd*odd[g(k,j-1)][j-1+m])*R;
	}
      }
  }  
  if (m%2 ==1) 
    return odd[0][m];
  else
    return even[0][m];
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

  double bars[20000];

  ofstream fp;
  fp.open("bfsg.dat");


  for( steps=10; steps<=50; steps+=10)
  {
    double dt = T/steps;
    double dw = sigma*sqrt(dt); // up movement
    double du = (r-sigma*sigma/2)*dt; // up per unit time
 
    int i,j, n;
    n=0;
    for(i=0; i<= steps;i++) {
      bars[n] = i*du+i*dw;
      n++;
      for(j=0;j<i;j++) {
	bars[n]=bars[n-1]-2*dw;
	n++;
      }
    }
    sort(bars,bars+n);

    double valp, valn;
    valp = 0;
    valn = exp(-r*T);
    value = 0;
    double barstart = log(X/S);
    for(i=0; i<n; i++) {
      if(bars[i]< barstart) continue;
      B= S*exp(bars[i]);
      valp = valn;
      valn=forwardST(S, r, B, alpha, sigma, T, steps); 
      //cerr<<i <<" "<< B<< " "<< valp << " " << valn << " " << value << endl;
      value=value+max(B-X,0.0)*(valp-valn);
      if(valp <= 0.00001) break;
    }
    cerr << valp << " " << i << endl;
    fp << steps << "  " << value << endl;
  }		
  fp.close();
  return 0;
}
