// Title: The Best / Rayner / Thas statistic for the Laplace distribution
// Ref. (book or article): Best, D. J., Rayner, J. C. W., & Thas, O. 2008. Comparison of some tests of fit for the Laplace
// distribution. Comput. Statist. Data Anal., 52(12), 5338–5343.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat104(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
      // Here, INDICATE the name of your statistic
      const char *nom = "$BRT_{3,4}$";
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
      double sqrt(double x);
      double statV3, statV4, xbar, s, b1sqrt, tmp = 0.0, tmp2 = 0.0, tmp3 = 0.0, tmp4 = 0.0, b2;
      
      // calculate mu^
      for (i = 0; i < n; i++) {
	tmp = tmp + x[i];
      }
      xbar = tmp / (double)n;
      
      // calculate b^
      for (i = 0; i < n; i++) {
	tmp2 = tmp2 + R_pow(x[i] - xbar, 2.0);
      }
      s = sqrt(tmp2 / (double)n);
            
      // calculate b1sqrt
      for (i = 0; i < n; i++) tmp3 = tmp3 + R_pow((x[i] - xbar) / s, 3.0);
      b1sqrt = tmp3 / (double)n;

      // calculate b2
      for (i = 0; i < n; i++) tmp4 = tmp4 + R_pow((x[i] - xbar) / s, 4.0);
      b2 = tmp4 / (double)n;
      
      // calculate statV3
      statV3 = b1sqrt * sqrt((double)n / 54.0);	
      statV4 = (b2 - 6.0) * sqrt((double)n / 1072.8);	

      statistic[0] = 6.0 * R_pow(statV3, 2.0) / 7.0 + 149.0 * R_pow(statV4, 2.0) / 165.0; // Here is the test statistic value	

      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue104.cpp"
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
