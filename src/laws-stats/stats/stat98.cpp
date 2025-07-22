// Title: The ANB statistic for the Laplace distribution
// Ref. (book or article): Tests of goodness of fit based on Phi-divergence
// Hadi Alizadeh Noughabi and Narayanaswamy Balakrishnan.
// Journal of Applied Statistics (2016)
// 43:3, 412-429

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat98(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$AB_{l}$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	if((int)paramstat[0] == 2) {
	  paramstat[0] = 2.0;
	  nom = "$AB_{He}$";
	} else if((int)paramstat[0] == 3) {
	  paramstat[0] = 3.0;
	  nom = "$AB_{Je}$";
	} else if((int)paramstat[0] == 4) {
	  paramstat[0] = 4.0;
	  nom = "$AB_{TV}$";
	} else if((int)paramstat[0] == 5) {
	  paramstat[0] = 5.0;
	  nom = "$AB_{\\chi^2}$";
	} else {
	  paramstat[0] = 1.0;
	  nom = "$AB_{KL}$";
	}
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
    int version;
    // version is: KL, He, Je, TV, Chi^2 (1 to 5)
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      version = 1;
      paramstat[0] = 1.0;
    } else if (nbparamstat[0] == 1) {
      version = (int)paramstat[0];
    } else {
      Rf_error("Number of parameters should be at most: 1");
    }

    // If necessary, we check if some parameter values are out of parameter space
    if ((version != 1) && (version != 2) && (version != 3) && (version != 4) && (version != 5)) {
      Rf_warning("version should be 1, 2, 3, 4 or 5 in stat98!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

    if (n > 3) {
// Computation of the value of the test statistic
      void R_rsort(double* x, int n);
      double hstar;
      double sqrtpi = sqrt(M_PI);
      double sqrt2pi = sqrt(2.0) * sqrtpi;
      
      R_rsort(x, n); 		// we sort the data

      double muhat, bhat, tmp;
      // calculate mu^ and b^ by using the maximum likelihood estimators 
      // mu^ = the sample median
      // b^ = 1/n * \sum_{i=1}^{n} |xi - mu^|
      
      // calculate mu^     
      if (n % 2 == 0) {		// check if n is divisible by 2
	muhat = (x[n / 2 - 1] + x[n / 2]) / 2.0;
      } else {
	muhat = x[n / 2];
      }
      
      // calculate b^
      tmp = 0.0;
      for (i = 0; i < n; i++) {
	tmp = tmp + fabs(x[i] - muhat);
      }
      bhat = tmp / (double)n;

      hstar = bhat * R_pow(2.0 / ((double)n * sqrtpi), 0.2);
      
      double zihat, D = 0.0, fhat, f, nu = 0.0;
      for (i = 1; i <= n; i++) {
	zihat = (x[i - 1] - muhat) / bhat;
	fhat = 0.0;
	for (j = 1; j <= n; j++) {
	  tmp = (x[j - 1] - x[i - 1]) / hstar;
	  fhat = fhat + exp(-0.5 * tmp * tmp);
	}
	fhat = fhat / ((double)n * hstar * sqrt2pi);
	f = 0.5 * exp(-fabs(zihat)) / bhat;
	tmp = fhat / f;
	if (version == 1) nu = tmp * log(tmp);
	if (version == 2) nu = 0.5 * R_pow(sqrt(tmp) - 1.0, 2.0);
	if (version == 3) nu = (tmp - 1.0) * log(tmp);
	if (version == 4) nu = fabs(tmp - 1.0);
	if (version == 5) nu = R_pow(tmp - 1.0, 2.0);
	D = D + nu / tmp;
      }
      D = D / (double)n;
      
      statistic[0] = D; // Here is the test statistic value
      

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue98.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i = 0; i <= (nblevel[0] - 1);i++) {
    if (usecrit[0] == 1) { // We use the provided critical values
	if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
		if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }
    
// If applicable, we free the unused array of pointers

}

// We return
    return;
   
        
  }
 
}
