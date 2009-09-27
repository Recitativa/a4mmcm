#include<iostream>
#include<list>
#include <cmath>
#include <algorithm>
#include <cassert>

using namespace std;

typedef double Real;

Real K; 
Real su, sdelta, dis, rho, alpha;
Real dt;


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
    OrderedPath.insert(it,Z);						\
    NTERM.in_it = (--it);						\
    assert(*(NTERM.in_it)==Z);						\
    /* Compute the quantile */						\
    it = PTERM.qu_it;							\
    int qth = PTERM.qth;						\
    if(Z<=*(PTERM.qu_it)) qth++;					\
    int di = qth - (int)(alpha*nstack);					\
    NTERM.qth = qth - di;						\
    if(di>0) 								\
      for(;di>0;di--,it--) NULL;					\
    if(di<0)								\
      for(;di<0;di++,it++) NULL;					\
    assert(NTERM.qth==(int)(alpha*nstack));				\
    NTERM.qu_it = it;							\
    Real Z1 = *it, Z2 = *(++it);					\
    /*   NTERM.Walpha = Z1 + (Z2-Z1)*(nstack*alpha-qth);*/		\
    /*NTERM.g1 = payfun(NTERM.Walpha);*/				\
    NTERM.Walpha = *it;							\
    if(nstack==steps-1) NTERM.g1 = payfun(NTERM.Walpha);       		\
    else NTERM.g1=1000000;						\
    nstack++;								\
    /*PRINT_LIST;*/							\
  } while(false)

  do {
    Node &PTERM = STACK[nstack-1]; // present Node;
    Node &NTERM = STACK[nstack]; // next Node;
    Real Z;
    //PRINT_LIST;
    if(nstack==steps)
      { pay = PTERM.g1;
	DELTERM;
      }
    else if(PTERM.st == down) {
      //pay = max(PTERM.g1, dis*(PTERM.g2+pay)/2);
      pay = dis*(PTERM.g2+pay)*.5;
      DELTERM;
    }
    else if(PTERM.st == up) {
      PTERM.g2 = pay;
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
  return pay;
}

typedef struct {
  Real h;
  Real A;
} RichNode;

Real richardson_extrapolation(RichNode * data, int len) {
  RichNode *S = new RichNode[len];
  int i,j,k;
  for(i=0;i<len;i++) S[i] = data[i];

  for(k=1,j=len;j>=2;j--,k++) {
    for(i=0;i<j-1;i++) {
      Real t = S[i].h/S[i+1].h;
      t = pow(t,k);
      S[i].A = (t*S[i+1].A-S[i].A)/(t-1);
      cerr << S[i].A << " ";
    }
    cerr << endl;
  }
  return S[0].A;
}


int main()
{
  Real S0, K, r, sigma,T, mu;
  int n,i;
  S0=100, K=95, alpha=0.8, r=0.05, sigma=0.2, T=.25, mu=r-sigma*sigma/2;
  n = 20;
  
  RichNode R[100];
  R[0].h = .025, R[0].A=10.2429;
  R[1].h = 0.0208333, R[1].A=9.79268;
  R[2].h = 0.0178571, R[2].A=9.48657;
  R[3].h = 0.015625, R[3].A=9.45113;
  R[4].h = 0.0138889, R[4].A=9.19926 ;
  cerr << richardson_extrapolation(R,5) << endl;;

  return 0;
  

  for(i=0;i<5;i++) {
    n= 10+ i*2;

    dt = T/n;
    su = mu*dt;// su=.5;
    sdelta= sigma*sqrt(dt); //sdelta=1;
    dis = exp(-r*dt);
    rho = S0/K;
    cout<< su << " " << sdelta <<" "<< alpha << " " << dis <<" " << rho << endl;
    OrderedPath.clear();
    R[i].h = dt;
    R[i].A = K*pprice(n+1,0);
    cout << dt <<" " << n << " " << R[i].A << endl;
  }
  cerr << "RICH" << endl;
  cout << richardson_extrapolation(R,i) << endl;
  return 0;
}



