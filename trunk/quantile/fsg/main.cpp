#include <iostream>
#include <fstream>
#include <cmath>     //standard mathematical library
#include <algorithm>
#include <vector>
#include <cassert>
using namespace std;

double factor=0.5;

double forward_shorting_grid_cumulative_parisian_binary_option( double S, //spot price
								double r, //interest rate
								double B, //down barrier
								double alpha, //
								double sigma,  //volatility	
								double T,  //time to maturity 
								int steps) {
  double dt = T/steps;
  double R=exp((-r)*(T/steps)); 	//discount rate at each step
  double up=sigma*sqrt(dt)/factor; // up movement
  double D=alpha*(steps+1);  //barrier
  double u=(sigma*sigma*dt + dt*dt*(r-sigma*sigma/2))/(up*up);
  double c=(r-sigma*sigma/2)*dt/up;
  double pu=(u+c)/2.;  //probability of upward
  double pd=(u-c)/2.;   //downward
  double p0=1.-u;      //zero
//  cout<<pu<<' '<<pd<<' '<<p0<<' '<<endl;
  int m=steps;
  double b=log(B/S)/up;
  double odd[m+1][2*m+1];
  double even[m+1][2*m+1];

#define g(k,j) (((j)<=b)?((k)+1):(k))
	
  for(int k=0; k<m+1;k++)
    for(int j=-m; j<m+1;j++) {
      if (k<D) even[k][j+m]=1;
      else even[k][j+m]=0;
    }
	
  for(int i=m-1;i >=0; i--)
   { for(int k=0; k<i+1;k++)
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
  vector<double> valp(2000*2+2,0);
  double value;
  
  ofstream fp;
  //fp.open("fsg.dat", ofstream::app);
  fp.open("fsg.dat");
  valp[0]=exp(-r*T);
  for( steps=100; steps<=150; steps+=10)
  {
    double up=sigma*sqrt(T/steps)/factor; // up movement
    for(int j= -steps; j<steps+1; j++) {
      B=S*exp(j*up);
      valp[j+steps+1] =							\
	forward_shorting_grid_cumulative_parisian_binary_option(S, r, B, alpha, sigma, T, steps); 
    }
 
    value=0;	
    for(int j=-steps; j<steps+1;j++)
      value=value+max(exp(j*up)*S-X,0.0)*(valp[j+steps]-valp[j+steps+1]);
    fp << steps << "  " << up << " " <<  value << endl;
    cerr << steps << " " << value << endl;
  }		
  fp.close();
  return 0;
}
