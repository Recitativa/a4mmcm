#include<iostream>
#include<list>
#include <cmath>
#include <algorithm>

using namespace std;

typedef double Real;

Real K; 
Real mu, sdelta, dis;
//Real (*payfun)();
Real alpha;

list<Real> OrderedPath;

Real payfun(Real Walpha) {
  return max(exp(Walpha)-1., .0);
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

Real pprice(int steps, Real ZZ) {
  list<Real>::iterator it;
  Real g1,g2,g3, pay;
  Real Walpha, Z;
  
  Node *STACK = new Node[steps+1];
  int nstack = 1;
  OrderedPath.push_back(Z);
  STACK[0].in_it = OrderedPath.begin();
  STACK[0].qu_it = OrderedPath.begin();
  STACK[0].st = no;
  STACK[0].qth = 0;
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
    }                                                                   \
    if(Z>*it) it++;							\
    OrderedPath.insert(it,Z);						\
    NTERM.in_it = it;							\
    it = PTERM.qu_it;							\
    int qth = PTERM.qth;						\
    if(Z<=*(PTERM.qu_it)) qth++;					\
    int di = qth - (int)(nstack*alpha);					\
    NTERM.qth = qth - di;						\
    if(di>0) 								\
      for(;di>0;di--) it--;						\
    if(di<0)								\
      for(;di<0;di++) it++;						\
    NTERM.qu_it = it;							\
    Real Z1 = *it, Z2 = *(++it);					\
    NTERM.Walpha = Z1 + (Z2-Z1)*(nstack*alpha-qth);			\
    NTERM.g1 = payfun(NTERM.Walpha);					\
    nstack++;								\
  } while(false)

  do {
    Node &PTERM = STACK[nstack-1]; // present Node;
    Node &NTERM = STACK[nstack]; // next Node;
    if(nstack==steps)
      { pay = PTERM.g1;
	DELTERM;
	nstack--;
      }
    else if(PTERM.st == down) {
      pay = max(PTERM.g1, dis*(PTERM.g2+pay)/2);
      DELTERM;
      nstack--;
    }
    else if(PTERM.st == up) {
      PTERM.g2 = pay;
      PTERM.st = down;
      Z = *(PTERM.in_it) + mu - sdelta;
      ADDTERM;
    }
    else if(PTERM.st == no){
      PTERM.st = up;
      Z = *(PTERM.in_it) + mu + sdelta;
      ADDTERM;
    } 
  }while(nstack==0);
  return pay;
}


int main()
{
  mu = 1;
  sdelta= 2;
  alpha = .5;
  dis = exp(-.01);
  // payfun=callfun; 
  OrderedPath.clear();
  cout << pprice(30,0) << endl; 
  
  return 0;
}



