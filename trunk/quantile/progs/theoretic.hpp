#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <cstdio>
#include <iostream>


using namespace std;

class Brownian {
public:
  double sigma;
  double mu;
  Brownian(double Isigma, double Imu): sigma(Isigma), mu(Imu) {};
};

typedef struct {
  double T; 
  int Terms; // How many loops times 
  int Rb; // records begin with 2^Rb+1 points 
  int Rm; // the quantile computation ends at 2^Rm
  int Re; // records end with 2^Re+1 points
  int Nseg;
  unsigned long int Rseed; // Random seed
} SimPara;


class BrownSim {
private:
  double sigma;
  double mu;
public:
  BrownSim(double Isigma = 1, double Imu=0): sigma(Isigma), mu(Imu) {}
  int Sim(SimPara Para);
};


class DifferentStepsException {};


// The step Iterator class
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
    return (p == rhs.p) && (step==rhs.step);}
  bool operator != (const StepIter<T> & rhs) const { return !((*this) == rhs); }
};
