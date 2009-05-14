#define BOOST_TEST_DYN_LINK
#define BOOST_AUTO_TEST_MAIN 
#include <iostream>
#include <boost/test/auto_unit_test.hpp> 
#include <cmath>

#include "theoretic.hpp"
#include <gsl/gsl_integration.h>

using namespace std;

double f(double x, void * params) {
  double alpha = *(double *) params;
  double ff = log(alpha*x) / sqrt(x);
  return ff;
}

double inte_test() {
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  double result, error;
  double expected = -4.0;
  double alpha = 1.0;
  
  gsl_function F;
  F.function = &f;
  F.params = &alpha;
  
  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
			w, &result, &error); 
  gsl_integration_workspace_free(w);
  return result;
}


BOOST_AUTO_TEST_CASE( test1 ) 
{ 
  // with or without some check 
  BOOST_REQUIRE( fabs(inte_test() + 4.0) < 0.01);
  
  Quatile Q1(1.0,.3);
  //cout << "g1(2,4,1,1)" << << endl;
  cout << "g(2,8) = " << Q1.g(1,2) << endl;

} 
