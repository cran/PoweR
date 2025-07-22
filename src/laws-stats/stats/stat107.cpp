// Title: Kozubowski / Panorska 
// Ref. (book or article): T.J. Kozubowski, A.K. Panorska / Journal of Statistical Planning and Inference 120 (2004) 41 -- 63

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat107(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
      // Here, INDICATE the name of your statistic
      const char *nom = "$KP$";
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
    
    if (n > 3) {
      // Computation of the value of the test statistic
      double Rf_pchisq(double q, double df, int lower_tail, int log_p);
      void R_rsort (double* x, int n);

      double muhat, kappahat, tmp, tmp1 = 0.0, tmp2 = 0.0;
      
      R_rsort(x, n); 		// we sort the data

      // mu^ = the sample median
      if (n % 2 == 0) {		// check if n is divisible by 2
	muhat = (x[n / 2 - 1] + x[n / 2]) / 2.0;
      } else {
	muhat = x[(n - 1) / 2];
      }

      for (i = 0; i < n; i++) {
	tmp = x[i] - muhat;
	if (tmp < 0) tmp1 = tmp1 - tmp; else tmp2 = tmp2 + tmp;
      }

      if (fabs(tmp2) < 0.0000000001) statistic[0] = 1.0; else {
	kappahat = tmp1 / tmp2;
	statistic[0] = (double)n * (2.0 - R_pow(1.0 + sqrt(kappahat), 2.0) / (1.0 + kappahat)); // Here is the test statistic value
      }
      
      
      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue107.cpp"
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
      
    }
    
    // We return
    return;
       
  }
  
}
