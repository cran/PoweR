// Title: The AN statistic for the Laplace distribution
// Ref. (book or article): Empirical Likelihood Ratio-Based Goodness-of-Fit Test
// for the Laplace Distribution
// Hadi Alizadeh Noughabi
// Commun. Math. Stat. (2016) 4:459–471


#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat99(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$A_{ratio}$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	paramstat[0] = 0.5;
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
    double delta;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      delta = 0.5;
      paramstat[0] = 0.5;
    } else if (nbparamstat[0] == 1) {
      delta = paramstat[0];
    } else {
      Rf_error("Number of parameters should be at most: 1");
    }

    // If necessary, we check if some parameter values are out of parameter space
    if ((delta <= 0.0) || (delta >= 1)) {
      Rf_warning("delta should be in (0, 1) in stat99!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }


    if (n > 3) {
// Computation of the value of the test statistic
      void R_rsort(double* x, int n);
      double term1, temp1, temp2;
      double *term2;
      term2 = new double[n];
      int m;
      double eps = 2.22044604925031308085 * R_pow(10.0, -16.0);
      double bound;
      int cte;
      
      R_rsort(x, n); 		// we sort the data

      double muhat, bhat, tmp = 0.0, min = x[0], max = x[n - 1];
      // calculate mu^ and b^ by using the maximum likelihood estimators 
      // mu^ = the sample median
      // b^ = 1/n * \sum_{i=1}^{n} |xi - mu^|
      
      // calculate mu^     
      if (n % 2 == 0) {		// check if n is divisible by 2
	muhat = (x[n / 2 - 1] + x[n / 2]) / 2.0;
      } else {
	muhat = x[(n - 1) / 2];
      }

      
      
      // calculate b^
      for (i = 0; i < n; i++) {
	tmp = tmp + fabs(x[i] - muhat);
      }
      bhat = tmp / (double)n;

      
      
      for (j = 1; j <= n; j++) {
	term2[j] = 0.5 * exp(-fabs(x[j - 1] - muhat) / bhat) / bhat;
      }


      tmp = 1.0; // m = 1 here
      for (j = 1; j <= n; j++) {
	if ((j + 1) <= n) temp1 = x[j]; else temp1 = max;
	if (j >= 2) temp2 = x[j - 2]; else temp2 = min;
	tmp = tmp * 2.0 / ((temp1 - temp2) * (double)n);
	tmp = tmp / term2[j];
      }

      term1 = tmp;

      bound = R_pow((double)n, delta);
      if (bound > ((double)n / 2.0)) bound = (double)n / 2.0;
      cte = (int)bound;
      if ((bound - (double)(int)bound) < eps) cte = cte - 1;
      for (m = 2; m <= cte; m++) {
	tmp = 1.0;
	for (j = 1; j <= n; j++) {
	  if ((j + m) <= n) temp1 = x[j + m - 1]; else temp1 = max;
	  if ((j - m) >= 1) temp2 = x[j - m - 1]; else temp2 = min;
	  tmp = tmp * (double)(2 * m) / ((temp1 - temp2) * (double)n);
	  tmp = tmp / term2[j];
	}
	if (tmp < term1) term1 = tmp;  
      }
	

      statistic[0] = term1; // Here is the test statistic value
      

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue99.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i = 0; i <= (nblevel[0] - 1); i++) {
    if (usecrit[0] == 1) { // We use the provided critical values
	if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
		if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }
    
// If applicable, we free the unused array of pointers
    delete[] term2;
    
}

// We return
    return;
   
        
  }
 
}
