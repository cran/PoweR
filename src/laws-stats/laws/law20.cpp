// Title: JSB(g,d)
// Ref. (book or article): N. L. Johnson, "Systems of Frequency Curves Generated by Methods of Translation", Biometrika, Vol. 36, No. 1/2 (Jun., 1949), pp. 149-176
//                         Kallenberg, Wilbert C. M. and Ledwina, Teresa (1997), "Data driven smooth tests for composite hypotheses: comparison of powers", 
//						   J. Statist. Comput. Simulation 59 (1997), No. 2, 101--121. 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void law20 (int *xlen, double *x, char **name, int *getname, double *params, int *nbparams, int *setseed) {

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of the distribution:
      const char *nom = "$JSB(g,d)$";
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
      for (i = j; i < 50; i++) name[i][0] = space[0];
      return;
    }
    
// Initialization of the parameters
    double g, d;
    if (nbparams[0] == 0) {
      nbparams[0] = 2;
      g = 0.0;
      d = 1.0;
      params[0] = 0.0;
      params[1] = 1.0;
    } else if (nbparams[0] == 1) {
      nbparams[0] = 2;
      g = params[0];
      d = 1.0;
      params[1] = 1.0;
    } else if (nbparams[0] == 2) {
      g = params[0];
      d = params[1];
    } else {
      Rf_error("Number of parameters should be at most: 2");
    }
   
// If necessary, we check if some parameter values are out of parameter space
    if (d <= 0.0) {
      Rf_warning("d should not be <= 0 in law20!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

// Generation of the random values
    if (setseed[0] == 1) GetRNGstate();
    double Rf_rnorm(double mu, double sigma);
    double *Z;
    Z = new double [n];
    for (i = 0; i < n; i++) Z[i] = Rf_rnorm(0.0, 1.0);
    for (i = 0; i < n; i++) x[i] = 1.0 / (1.0 + exp(-(Z[i] - g) / d));
    if (setseed[0] == 1) PutRNGstate();

// If applicable, we free the unused array of pointers. Then we return.
    delete[] Z;
    return;
    
  }
  
}
