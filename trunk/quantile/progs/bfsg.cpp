#include <iostream>
#include <fstream>
#include <sstream>
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

int main(int argc,      
          char *argv[],  
          char *envp[] )
{
  double S; //spot price
  double K;  //strike price
  double r; //interest rate
  double alpha; 
  double sigma;  //volatility
  double T;  //time to maturity 

  int Bn, En;
 
  if(argc != 9) {cerr << argv[0] <<" S K alpha r sigma T Bn En" << endl; return 1;}
  stringstream ss (stringstream::in | stringstream::out);
  for(int i=1;i<9;i++) { ss << argv[i] << " ";}
  ss >> S >> K >> alpha >> r >> sigma >> T >> Bn >> En;

  ostringstream SoutFilename;
  SoutFilename << "bfsg_S" << S << "_K"<< K << "_alpha"<< alpha << "_r" << r << "_sigma_" << sigma << "_T_" << T<< "_Bn_"<< Bn <<"_En_"<< En;
  SoutFilename << ".txt";

  ofstream outf;
  
  outf.open(SoutFilename.str().c_str(),ios::app);
  outf << "Pricing S0:" << S << " K:"<< K << " alpha:"<< alpha << " r:" << r << " sigma:" << sigma << " T:" << T << endl;


  double B; //down barrier
  int steps;
  double value;

  double *bars=new double[1024*1024];

  for( steps=Bn; steps<=En; steps+=10)
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
    double barstart = log(K/S);
    for(i=0; i<n; i++) {
      if(bars[i]< barstart) continue;
      B= S*exp(bars[i]);
      valp = valn;
      valn=forwardST(S, r, B, alpha, sigma, T, steps); 
      //cerr<<i <<" "<< B<< " "<< valp << " " << valn << " " << value << endl;
      value=value+max(B-K,0.0)*(valp-valn);
      if(valp <= 0.00001) break;
    }
    cerr << valp << " " << i << endl;
    outf << steps << "  " << value << endl;
  }		
  outf.close();
  return 0;
}
