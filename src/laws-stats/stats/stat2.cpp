// Title: The Anderson-Darling test
// Ref. (book or article): See package nortest and also Table 4.9 p. 127 in M. A. Stephens, “Tests Based on EDF Statistics,” In: R. B. D’Agostino and M. A. Stephens, Eds., Goodness-of-Fit Techniques, Marcel Dekker, New York, 1986, pp. 97-193.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat2(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$AD^*$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 0;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	
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
    
    if (n > 7) {
      // Computation of the value of the test statistic
      //    double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
      void R_rsort (double* x, int n);
      double *logp1, *logp2;
      logp1 = new double[n];
      logp2 = new double[n];
      double statAD, varX, sdX, meanX, sumAD, pval, tmp;
      R_rsort (x, n); // We sort the data
      
      // Two-pass algorithm algorithm is less prone to numerical errors than the naive approach
      // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
      meanX = 0.0;
      for (i = 0; i < n; i++) meanX = meanX + x[i];
      meanX = meanX / (double)n;
      varX = 0.0;
      for (i = 0; i < n; i++) varX = varX + R_pow(x[i] - meanX, 2.0);
      varX = varX / (double)(n - 1); 
      sdX = sqrt(varX);
      
      for (i = 0; i < n; i++) {
	tmp = (x[i] - meanX) / sdX;
        logp1[i] = Rf_pnorm5(tmp, 0.0, 1.0, 1, 1);
        logp2[i] = Rf_pnorm5(-tmp, 0.0, 1.0, 1, 1);
      }
    
      sumAD = 0.0;
      for (i = 1; i <= n; i++) sumAD = sumAD + (double)(2 * i - 1) * (logp1[i - 1] + logp2[n - i]);
      statAD = -((double)n + sumAD / (double)n);
    
      statistic[0] = statAD; // Here is the test statistic value
      
      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue2.cpp"
      }
      
// We take the decision to reject or not to reject the null hypothesis H0
      for (i = 0; i < nblevel[0]; i++) {
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
	} else {
	  if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
	}
      }
    
// If applicable, we free the unused array of pointers
      delete[] logp1;
      delete[] logp2;

}

// We return
    return;
   
        
  }
  
}
