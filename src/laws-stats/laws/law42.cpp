// Title: Inverse Gaussian(mu,lambda)
// Ref. (book or article): https://cran.r-project.org/package=statmod
// Giner, G., and Smyth, G. K. (2016). statmod: Probability calculations for the inverse Gaussian
// distribution. R Journal 8(1), 339-351. https://journal.r-project.org/archive/2016-1/giner-smyth.pdf

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law42(int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$invgauss(mu,lambda)$";
// Here, INDICATE the number of parameters of the distribution:
      nbparams[0] = 2;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe params has not be initialized with a sufficient length since the correct value of nbparams[0] may be unkown yet).
	params[0] = 1.0;
	params[1] = 1.0;
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
    double mu, lambda;
    if (nbparams[0] == 0) {
      nbparams[0] = 2;
      mu = 1.0;
      lambda = 1.0;
      params[0] = 1.0;
      params[1] = 1.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 2;
      mu = params[0];
      lambda = 1.0;
      params[1] = 1.0;
    } else if  (nbparams[0] == 2) {
      mu = params[0];
      lambda = params[1];
   } else {
      Rf_error("Number of parameters should be at most: 2");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (mu <= 0.0 || lambda <= 0.0) {
      Rf_warning("mu and lambda should be > 0 in law42!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();
    double rf_rnorm(double mu, double sigma);
    double Rf_runif(double a, double b);
    double bignumber = 5.0 * R_pow(10.0, 5.0);
    double x1;
    double phi = 1.0 / lambda;
 
    double *Y;
    Y = new double [n];
    // I used N(0,1)^2 instead of Chi^2(1) here to get exam same results as statmod::rinvgauss
    for (i = 0; i < n; i++) Y[i] = R_pow(Rf_rnorm(0.0, 1.0), 2.0) * phi * mu;

    for (i = 0; i < n; i++)   {
      if (Y[i] > bignumber) {
	x1 = 1.0 / Y[i];
	if (Rf_runif(0.0, 1.0) < 1 / (1 + x1)) x[i] = mu * x1; else x[i] = mu / x1;
      } else {
	x1 = 1.0 + 0.5 * Y[i] * (1.0 - sqrt(1.0 + 4.0 / Y[i]));
	if (Rf_runif(0.0, 1.0) < 1 / (1 + x1)) x[i] = mu * x1; else x[i] = mu / x1;
      }
    }

    
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] Y;

  return;
    
  }
  
}
