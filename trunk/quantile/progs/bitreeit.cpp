#include <omp.h>
#include <iostream>
#include <list>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>
#include "cdf.hpp"

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
  Real pprice(int steps, Real ZZ);
};

#define PRINT_LIST       do {				\
    list<Real>::iterator iit;				\
    for(iit=OrderedPath.begin();			\
	iit!=OrderedPath.end();				\
	iit++)						\
      cout << *iit << " ";				\
    cout << endl;					\
  }while(false)					



enum node_stat {no, up, down};

typedef struct{
  list<Real>::iterator in_it;
  list<Real>::iterator qu_it;
  node_stat st;
  Real pay, vup,vdown;
  Real Walpha;
  int qth;
} Node;


Real Option::pprice(int steps, Real ZZ) {
  list<Real>::iterator it;
  Real pay;
  Real Walpha;
  
  Node *STACK = new Node[steps+1];
  int nstack = 1;

  OrderedPath.clear();
  OrderedPath.push_front(ZZ);
  STACK[0].in_it = OrderedPath.begin();
  STACK[0].qu_it = OrderedPath.begin();
  STACK[0].st = no;
  STACK[0].qth = 0;

  //PRINT_LIST;

#define DELTERM do {				\
    OrderedPath.erase(PTERM.in_it);		\
    nstack--;					\
  } while(false)
  
#define ADDTERM do {							\
    NTERM.st = no;							\
    /* Insert node to OrderedPath*/					\
    if(Z>*(PTERM.in_it)) {						\
      for(it=PTERM.in_it;it!=OrderedPath.end() && Z > (*it);it++)	\
	NULL;								\
    }									\
    else {								\
      for(it=PTERM.in_it;it!=OrderedPath.begin() && Z <=*(it);it--)	\
	NULL;								\
      if(Z>*it) it++;							\
    }									\
    NTERM.in_it = OrderedPath.insert(it,Z);				\
    assert(*(NTERM.in_it)==Z);						\
    /* Compute the quantile */						\
    it = PTERM.qu_it;							\
    int qq = PTERM.qth;						\
    if(Z<=*(PTERM.qu_it)) qq++;					\
    int di = qq - (int)(alpha*nstack);					\
    NTERM.qth = qq - di;						\
    for(;di>0;di--,it--) NULL;						\
    for(;di<0;di++,it++) NULL;						\
    assert(NTERM.qth==(int)(alpha*nstack));				\
    NTERM.qu_it = it;							\
    Real Z1 = *it, Z2 = *(++it);					\
    NTERM.Walpha = Z1 + (Z2-Z1)*(alpha*nstack-NTERM.qth);		\
    /*NTERM.Walpha = *it;*/						\
    NTERM.pay = payfun(NTERM.Walpha);					\
    /*PRINT_LIST;*/							\
    nstack++;								\
  } while(false)

  do {
    Node &PTERM = STACK[nstack-1]; // present Node;
    Node &NTERM = STACK[nstack]; // next Node;
    Real Z;
    //PRINT_LIST;
    if(nstack==steps+1)
      { pay = PTERM.pay;
	DELTERM;
      }
    else if(PTERM.st == down) {
      //pay = max(PTERM.g1, dis*(PTERM.g2+pay)*.5);
      pay = dis*(PTERM.vup+pay)*.5;
      DELTERM;
    }
    else if(PTERM.st == up) {
      PTERM.vup = pay;
      PTERM.st = down;
      Z = *(PTERM.in_it) + su - sdelta;
      ADDTERM;
    }
    else if(PTERM.st == no) {
      PTERM.st = up;
      Z = *(PTERM.in_it) + su + sdelta;
      ADDTERM;
    } 
  }while(nstack!=0);
  assert(OrderedPath.empty());
  return pay;
}

Real Option::EPrice(Real S0, Real K, Real T, int n)
{
  mu=r-sigma*sigma/2;
  dt = T/n;
  su = mu*dt;
  sdelta= sigma*sqrt(dt);
  dis = exp(-r*dt);
  rho = S0/K;
  OrderedPath.clear();
  return K*pprice(n,0);
}

int main(int argc,      
          char *argv[],  
          char *envp[] )
{
  Real S0, K, alpha, r, sigma, T;
  int n,i;
  int Bn, En;
 
  if(argc != 9) {cerr << argv[0] <<" S0 K alpha r sigma T Bn En" << endl; return 1;}
  stringstream ss (stringstream::in | stringstream::out);
  for(i=1;i<9;i++) { ss << argv[i] << " ";}
  ss >> S0 >> K >> alpha >> r >> sigma >> T >> Bn >> En;

  //  S0=100, K=95, alpha=0.5, r=0.05, sigma=0.2, T=1, 
  //Bn=20;En=40;


  ostringstream SoutFilename;
  SoutFilename << "biit_S0" << S0 << "_K_"<< K << "_alpha"<< alpha << "_r" << r << "_sigma_" << sigma << "_T_" << T<< "_Bn_"<< Bn <<"_En_"<< En;
  SoutFilename << ".txt";

  
  ofstream outf;
  
  outf.open(SoutFilename.str().c_str(),ios::app);
  outf << "Pricing S0:" << S0 << " K:"<< K << " alpha:"<< alpha << " r:" << r << " sigma:" << sigma << " T:" << T << endl;

#pragma omp parallel shared(outf)  num_threads(10)
  {
#pragma omp single
    {
      cerr << "num threads:"<<omp_get_num_threads()<<endl;
    }
#pragma omp for
    for(n=Bn;n<= En; n+=1)
      {
	Option A(r, sigma, alpha);
	Real price;
	price = A.EPrice(S0,K,T,n);
	outf << n << " " << price << endl;  
      }
  }

  outf.close();
  return 0;
}



