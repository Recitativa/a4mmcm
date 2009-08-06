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
list<Real>::iterator tit1,tit2;


Real payfun(Real Walpha) {
  return max(exp(Walpha)-1., .0);
}


Real pprice(int steps, Real Z) {
  list<Real>::iterator it;
  Real g1,g2,g3, pay;
  
  for(it=OrderedPath.begin();it!=OrderedPath.end() && Z< (*it);it++)
    NULL;
  OrderedPath.insert(it,Z);
  --it;

  Real Walpha;
  Real qq = (OrderedPath.size()*alpha);
  int nq = (int)qq;
  int i;
  for(tit1=OrderedPath.begin(),i=0;
      tit1!=OrderedPath.end() && i<nq; 
      tit1++, i++) 
    NULL;
  Real Z1 = *tit1, Z2 = ++tit1!=OrderedPath.end()?*(tit1):Z1;
  Walpha = Z1 + (Z2-Z1)*(qq-nq);


  g1 = payfun(Walpha);
  if(steps==0) return g1;
  g2 = pprice(steps-1, Z+mu+sdelta);
  g3 = pprice(steps-1, Z+mu-sdelta);
  pay = max(g1, dis*(g2+g3)/2);
  OrderedPath.erase(it);
  return pay;
};

int main()
{
  mu = 0.4;
  sdelta=.2;
  alpha = .5;
  dis = exp(-.01);
  // payfun=callfun; 
  OrderedPath.clear();
  cout << pprice(10,0) << endl; 
  
  return 0;
}



