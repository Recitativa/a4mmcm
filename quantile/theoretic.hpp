
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <cstdio>
#include <iostream>

//#include <stl_relops.h>

using namespace std;

// path-dependent options, The annals of Applied Probability, 1995,
// Vol. 5, No. 2, 389-298

const double PI = 3.14159265358979323846;

class Brownian {

public:
  double sigma;
  double mu;
  Brownian(double Isigma, double Imu): sigma(Isigma), mu(Imu) {};
};

class Quantile : public Brownian {

public:
  double alpha;
  Quantile(double Isigma = 1, double Ialpha = 1)
      : Brownian(Isigma, 0), alpha(Ialpha) {}

  double dX(double x);
  double mean();

private:
  static double dX_i(double x, void * Iparams);
  static double XdX_i(double x, void * Iparams) {
    return x*dX_i(x, Iparams);
  }
};



// Quantile Range Counter

class OutofRangeException {
public:
  size_t i;
  OutofRangeException(size_t I): i(I) {}
};

template<class Num, class C>
class QRCounter {
private:
  C total;
  Num Rbegin; //Begin of the range
  Num Rend; // End of the range
  Num h; // Step length
  size_t n; // number of point between Rbegin and Rend.
  C * counter;
public:
  QRCounter(Num begin, Num end, size_t In);
  ~QRCounter();
  int Add(Num x);
  void PrintDest();
  Num QuantileC(double alpha);
  Num Quantile(double alpha);
  int init();
};

template<class Num, class C>
QRCounter<Num, C>::QRCounter(Num begin, Num end, size_t In):
  total(0), Rbegin(begin), Rend(end), h((end-begin)/In),n(In) {
  counter = new C[n];
  init();
}
template<class Num, class C> 
int QRCounter<Num, C>::init() {
  total = 0;
  memset(counter, 0, sizeof(C)*n);
  return 0;
}

template<class Num, class C>
QRCounter<Num, C>::~QRCounter() {
  delete counter;
}

template<class Num, class C>
int QRCounter<Num, C>::Add(Num x) {
  long i = (size_t)((x-Rbegin)/h);
  if( i <0 ) { i=0;}
  if(i > (long)(n-1)) {i=(long)n-1;}
  counter[(size_t)i]++;
  total++;
  return 0;
}

template<class Num, class C>
void QRCounter<Num, C>::PrintDest() {
  std::cout << "n:" << n << std::endl;
  for(size_t i =0; i< n; i++) {
    std::cout << "[" << h*i+Rbegin << ":" << h*(i+1)+Rbegin << "]= "<<  counter[i] << " ";
  }
}

// Countinous Quantile
template<class Num, class C>
Num QRCounter<Num, C>::QuantileC(double alpha) {
  size_t i,j;
  C aim = (C)(alpha* total);
  C r1 = 0;
  C r2 = 0;
  for(i=0,j = 0; i< n; i++) {
    r1 = r2;
    r2 += counter[i];
    if (r2 > aim) break;
    if (counter[i]>0) j=i;
  }
  if (i==0 || i == n-1) 
    throw OutofRangeException(i);

  Num result = Rbegin+h*(Num)(j) + h*(Num)(aim-r1) * (Num)(j-i) /(Num)(r2-r1);
  return result;
}

// Quantile = inf {x #{X(i) < x} > alpha*T}
template<class Num, class C>
Num QRCounter<Num, C>::Quantile(double alpha) {
  size_t i;
  C aim = (C)(alpha* total);
  C r1 = 0;
  C r2 = 0;
  for(i=0; i< n; i++) {
    r1 = r2;
    r2 += counter[i];
    if (r2 > aim) break;
  }
  if (i==0 || i == n-1) 
    throw OutofRangeException(i);

  Num result = Rbegin+h*(Num)(i+1);
  return result;
}

class BrownSim {
private:
  double sigma;
  double mu;
public:
  BrownSim(double Isigma = 1, double Imu=0): sigma(Isigma), mu(Imu) {}
  void Sim(double, int);
};


class DifferentStepsException {};


#include <iterator>

template < class T>
struct StepIter :
  public iterator<std::random_access_iterator_tag, T, int> {
  private:
  T * p;
  int step;
  public:
  explicit StepIter(T * Ip=NULL, int Istep=0): p(Ip), step(Istep) {};
  StepIter<T>& operator= (T *Ip) {p = Ip; return *this;}
  StepIter<T> operator+ (int n) {
    StepIter<T> tmp=*this; tmp.p +=  n*step; return tmp;}
  StepIter<T>& operator+= (int n)  {p += n*step; return *this;}
  StepIter<T> operator- (int n) {
    StepIter<T> tmp=*this; tmp.p -=  n*step; return tmp;}
  StepIter<T>& operator-= (int n) {p -= n*step; return *this;}
  int operator- (const StepIter<T> & rhs ) {
    if(step==rhs.step) return (p-rhs.p)/step;
    else throw DifferentStepsException(); 
  }
  T& operator[] (int n) {return p[n*step];}
  T& operator* () {return *p;}
  StepIter<T>& operator--() {p -= step; return *this; }
  StepIter<T>& operator++() {p += step; return *this; }
  StepIter<T>  operator++(int) {
    StepIter<T> tmp = *this; ++(*this); return tmp; }
  StepIter<T>  operator--(int) {
    StepIter<T> tmp = *this; --(*this); return tmp; }
  bool operator< (const StepIter<T> & rhs) const { return (p-rhs.p)/step < 0;}
  bool operator > (const StepIter<T> & rhs) const { return (rhs< *this);}
  bool operator <= (const StepIter<T> & rhs) const { return !(rhs < *this);}
  bool operator >= (const StepIter<T> & rhs) const { return rhs <= *this;}
  bool operator== (const StepIter<T> & rhs) const {
    return p == rhs.p && step==rhs.step;}
  bool operator != (const StepIter<T> & rhs) const { return !(*this == rhs); }
};
