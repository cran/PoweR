// Title: Gumbel(mu,sigma)
// Ref. (book or article): Coles, Stuart(2001), "An introduction to statistical modeling of extreme values"
//	                       Springer Series in Statistics. Springer-Verlag London, Ltd., London, 2001. xiv+208 pp. ISBN: 1-85233-459-262-01 See page 48
//                         or see function rgumbel() in package VGAM 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law26 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$Gumbel(\\mu,\\sigma)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 2;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
      params[0] = 0.0;
      params[1] = 1.0;
     }
// The following 7 lines should NOT be modified
      const char *space = " ";
      while (nom[j] != '\0') {
	name[j][0] = nom[j];
	j++;
      }
      for (i = j;i < 50; i++) name[i][0] = space[0];
      return;
    }

// Initialization of the parameters
    double mu, sigma;
    if (nbparams[0] == 0) {
      nbparams[0] = 2;
      mu = 0.0;
      sigma = 1.0;
      params[0] = 0.0;
      params[1] = 1.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 2;
      mu = params[0];
      sigma = 1.0;
      params[1] = 1.0;
    } else if  (nbparams[0] == 2) {
      mu = params[0];
      sigma = params[1];
    } else {
      Rf_error("Number of parameters should be at most: 2");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (sigma < 0.0) {
      Rf_warning("sigma should not be < 0 in law26!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();
    double Rf_runif(double a, double b);
    for (i = 0; i < n; i++) x[i] = mu - sigma * log(-log(Rf_runif(0.0, 1.0)));
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    return;
 }

}



