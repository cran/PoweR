// Title: Skewt distribution(xi,omega,alpha,nu)
// Ref. (book or article): https://cran.r-project.org/package=sn

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law41(int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$Skewt(xi,omega,alpha,nu)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 4;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 0.0;
      params[3] = R_PosInf;
     }
// The following 7 lines should NOT be modified
      const char *space = " ";
      while (nom[j] != '\0') {
		name[j][0] = nom[j];
		j++;
      }
      for (i = j; i < 50; i++) name[i][0] = space[0];
      return;
    }

// Initialization of the parameters
    double xi, omega, alpha, nu;
    if (nbparams[0] == 0) {
      nbparams[0] = 4;
      xi = 0.0;
      omega = 1.0;
      alpha = 0.0;
      nu = R_PosInf;
      params[0] = 0.0;
      params[1] = 1.0;
      params[2] = 0.0;
      params[3] = R_PosInf;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 4;
      xi = params[0];
      omega = 1.0;
      alpha = 0.0;
      nu = R_PosInf;
      params[1] = 1.0;
      params[2] = 0.0;
      params[3] = R_PosInf;
    } else if  (nbparams[0] == 2) {
      nbparams[0] = 4;
      xi = params[0];
      omega = params[1];
      alpha = 0.0;
      nu = R_PosInf;
      params[2] = 0.0;
      params[3] = R_PosInf;
    } else if  (nbparams[0] == 3) {
      nbparams[0] = 4;
      xi = params[0];
      omega = params[1];
      alpha = params[2];
      nu = R_PosInf;
      params[3] = R_PosInf;
    } else if  (nbparams[0] == 4) {
      xi = params[0];
      omega = params[1];
      alpha = params[2];
      nu = params[3];
    } else {
      Rf_error("Number of parameters should be at most: 4");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (omega <= 0.0 || nu <= 0.0) {
      Rf_warning("omega and nu should be > 0 in law41!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();
    double Rf_rchisq(double df);
    double Rf_rnorm(double a, double b);
    double Rf_qnorm5(double p, double mean, double sd, int lower_tail, int log_p);
    double *u1, *u2, v;
    u1 = new double [n];
    u2 = new double [n];

    for (i = 0; i < n; i++) u1[i] = Rf_rnorm(0.0, 1.0);
    for (i = 0; i < n; i++) u2[i] = Rf_rnorm(0.0, 1.0);
    for (i = 0; i < n; i++) {
      if (u2[i] > alpha * u1[i]) {
	u1[i] = - omega * u1[i];
      } else {
	u1[i] = omega * u1[i];
      }
    }
  
    if(nu < R_PosInf) {  
      for (i = 0; i < n; i++) {
	v = Rf_rchisq(nu) / nu;
	x[i] = u1[i] / sqrt(v) + xi;
      }
    } else {
      for (i = 0; i < n; i++) x[i] = u1[i] + xi;
    }

    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] u1;
    delete[] u2;

  return;
    
  }
  
}
