#include <omp.h>
#include <iostream>
#include <list>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>


using namespace std;

typedef double Real;

class Option {
public:
  Option(Real rr, Real ssigma, Real aalpha ):
    r(rr),sigma(ssigma),alpha(aalpha) {};
  Real r,sigma,alpha; 
  ~Option() {
    OrderedPath.clear();
  }
  Real EPrice(Real S0, Real K, Real T, int n);
  Real payfun(Real Walpha) {
    return max(rho*exp(Walpha)-1., .0); }
  Real su, sdelta, dis, rho, mu;
  Real dt;
private:
  list<Real> OrderedPath;
  Real ppprice(int steps, Real ZZ, int n);
};


Real Option::ppprice(int steps, Real ZZ,int n) {
  Real g1,g2,Z;
  list<Real>::iterator it;
  OrderedPath.push_front(ZZ);
  it = OrderedPath.begin();
  if(steps==0) {
    Real Walpha,Z1,Z2;
    list<Real>::iterator iit;
    int i;
    OrderedPath.sort();
    for(iit = OrderedPath.begin(), i=0;
	i<n*alpha; i++,iit++) NULL;
    Z2 = *iit;
    Z1 = *(--iit);
    Walpha = Z1+(Z2-Z1)*(alpha*n-i+1);
    OrderedPath.erase(it);
    return payfun(Walpha);
  }
  g1 = ppprice(steps-1,ZZ+su+sdelta, n);
  g2 = ppprice(steps-1,ZZ+su-sdelta, n);
  OrderedPath.erase(it);
  return dis*(g1+g2)*.5;
}

Real Option::EPrice(Real S0, Real K, Real T, int n)
{
  Real p;
  mu=r-sigma*sigma/2;
  dt = T/n;
  su = mu*dt;
  sdelta= sigma*sqrt(dt);
  dis = exp(-r*dt);
  rho = S0/K;
  OrderedPath.clear();
  p= K*ppprice(n,0.0,n);
  assert(OrderedPath.empty());
  return p;
}


int main(int argc,      
          char *argv[],  
          char *envp[] )
{
  Real S0, K, alpha, r, sigma, T;
  int n,i;
  int Bn, En;
 
  if(argc != 9) {cerr << "bitree S0 K alpha r sigma T Bn En" << endl; return 1;}
  stringstream ss (stringstream::in | stringstream::out);
  for(i=1;i<9;i++) { ss << argv[i] << " ";}
  ss >> S0 >> K >> alpha >> r >> sigma >> T >> Bn >> En;

  //S0=100, K=95, alpha=0.5, r=0.05, sigma=0.2, T=1, 
  //Bn=20;En=40;

  ostringstream SoutFilename;
  SoutFilename << "bi_S0" << S0 << "_K_"<< K << "_alpha"<< alpha << "_r" << r << "_sigma_" << sigma << "_T_" << T<< "_Bn_"<< Bn <<"_En_"<< En;
  SoutFilename << ".txt";
  
  ofstream outf;
  
  outf.open(SoutFilename.str().c_str(),ios::app);
  outf << "Pricing S0:" << S0 << " K:"<< K << " alpha:"<< alpha << " r:" << r << " sigma:" << sigma << " T:" << T << endl;

  Option A(r, sigma, alpha);
 
  Real price;

#pragma omp parallel shared(outf)  num_threads(10)
  {
#pragma omp for
    for(n=Bn;n<= En; n+=1) {
      price = A.EPrice(S0,K,T,n);
      outf << n << " " << price << endl;  
    }
  }

  outf.close();
  return 0;
}
