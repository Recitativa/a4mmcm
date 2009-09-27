#include<iostream>
#include<list>
#include <cmath>
#include <algorithm>

using namespace std;

typedef double Real;

Real K; 
Real su, sdelta, dis, rho, alpha;
//Real (*payfun)();


list<Real> OrderedPath;

Real payfun(Real Walpha) {
  return max(rho*exp(Walpha)-1., .0);
}

enum node_stat {no, up, down};

typedef struct{
  list<Real>::iterator in_it;
  list<Real>::iterator qu_it;
  node_stat st;
  Real g1, g2;
  Real Walpha;
  int qth;
} Node;

#define PRINT_LIST       do {				\
    list<Real>::iterator iit;				\
    for(iit=OrderedPath.begin();			\
	iit!=OrderedPath.end();				\
	iit++)						\
      cout << *iit << " ";				\
    cout << endl;					\
  }while(false)					



Real pprice(int steps, Real ZZ) {
  list<Real>::iterator it;
  Real g1,g2,g3, pay;
  Real Walpha;
  
  Node *STACK = new Node[steps+1];
  int nstack = 1;
  cout << alpha << endl;
  OrderedPath.clear();
  OrderedPath.push_back(ZZ);
  STACK[0].in_it = OrderedPath.begin();
  STACK[0].qu_it = OrderedPath.begin();
  STACK[0].st = no;
  STACK[0].qth = 0;

  //PRINT_LIST;

#define DELTERM do {				\
    OrderedPath.erase(PTERM.in_it);		\
  } while(false)
  
#define ADDTERM do {							\
    NTERM.st = no;							\
    if(Z>*(PTERM.in_it)) {						\
      for(it=PTERM.in_it;it!=OrderedPath.end() && Z > (*it);it++)	\
	NULL;								\
    }									\
    else {						\
      for(it=PTERM.in_it;it!=OrderedPath.begin() && Z <=*(--it);)	\
	NULL;								\
      if(Z>*it) it++;							\
    }									\
    OrderedPath.insert(it,Z);						\
    NTERM.in_it = (--it);						\
    it = PTERM.qu_it;							\
    int qth = PTERM.qth;						\
    if(Z<=*(PTERM.qu_it)) qth++;					\
    int di = qth - (int)(alpha*nstack);					\
    NTERM.qth = qth - di;						\
    if(di>0) 								\
      for(;di>0;di--,it--) NULL;					\
    if(di<0)								\
      for(;di<0;di++,it++) NULL;					\
    NTERM.qu_it = it;							\
    Real Z1 = *it, Z2 = *(++it);					\
    NTERM.Walpha = Z1 + (Z2-Z1)*(nstack*alpha-qth);			\
    /*NTERM.g1 = payfun(NTERM.Walpha);*/				\
    if(nstack==steps-1) NTERM.g1 = payfun(NTERM.Walpha);       		\
    else NTERM.g1=10000;						\
    nstack++;								\
    /*PRINT_LIST;*/							\
  } while(false)

  do {
    Node &PTERM = STACK[nstack-1]; // present Node;
    Node &NTERM = STACK[nstack]; // next Node;
    Real Z;
    if(nstack==steps)
      { pay = PTERM.g1;
	DELTERM;
	nstack--;
      }
    else if(PTERM.st == down) {
      //pay = max(PTERM.g1, dis*(PTERM.g2+pay)/2);
      pay = dis*(PTERM.g2+pay)/2;
      DELTERM;
      nstack--;
    }
    else if(PTERM.st == up) {
      PTERM.g2 = pay;
      PTERM.st = down;
      Z = *(PTERM.in_it) + su - sdelta;
      ADDTERM;
    }
    else if(PTERM.st == no){
      PTERM.st = up;
      Z = *(PTERM.in_it) + su + sdelta;
      ADDTERM;
    } 
  }while(nstack!=0);
  return pay;
}


int main()
{

  Real S0, K, r, sigma,T, mu;
  int n;
  S0=100, K=95, alpha=0.8, r=0.05, sigma=0.2, T=.25, mu=r-sigma*sigma/2;
  n = 20;

  su = mu*T/n;
  sdelta= sigma*sqrt(T/n);
  dis = exp(-r*T/n);
  rho = S0/K;

  cout<< su << " " << sdelta <<" "<< alpha << " " << dis <<" " << rho << endl; 
  // payfun=callfun; 
  OrderedPath.clear();
  cout << K*pprice(n+1,0) << endl; 
  return 0;
}



