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
  
  Quantile Q1(1.0,.3);

  double dX1_E = 0.0612817;
  double mean_E =  -0.230539;
  double dX1 = Q1.dX(1);
  double mean = Q1.mean();
  cout << "dX(1) = " << dX1_E << endl;
  cout << "mean " << mean_E << endl;
  
  BOOST_REQUIRE( fabs(mean_E-mean) < 0.001);
  BOOST_REQUIRE( fabs(dX1_E -dX1) < 0.001);
  

} 
