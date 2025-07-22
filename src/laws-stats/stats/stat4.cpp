// Title: The Zhang-Wu Z_A test
// Ref. (book or article): Zhang, Jin and Wu, Yuehua (2005), Likelihood-ratio tests for normality,
//	                       Computational Statistics & Data Analysis, Vol. 49, No. 3, pp. 709-721. 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat4(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$Z_A$";
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
    void R_rsort (double* x, int n);
    //    double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
    double logp1, logp2;
    double varX, meanX, sdX, statZA, tmp;

    // Two-pass algorithm algorithm is less prone to numerical errors than the naive approach
    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    meanX = 0.0;
    for (i = 0; i < n; i++) meanX = meanX + x[i];
    meanX = meanX / (double)n;
    varX = 0.0;
    for (i = 0; i < n; i++) varX = varX + R_pow(x[i] - meanX, 2.0);
    varX = varX / (double)(n - 1); 
    sdX = sqrt(varX);

    R_rsort (x, n); // We sort the data

    statZA = 0.0;
    for (i = 1; i <= n; i++) {
      tmp = (x[i - 1] - meanX) / sdX;
      logp1 = Rf_pnorm5(tmp, 0.0, 1.0, 1, 1);
      logp2 = Rf_pnorm5(tmp, 0.0, 1.0, 0, 1);
      statZA = statZA - logp1 / ((double)(n - i) + 0.5) - logp2 / ((double)i - 0.5);
    }
    
    statZA = 10.0 * statZA - 32.0; // See p. 711 of the paper
    
    statistic[0] = statZA; // Here is the test statistic value


    if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue4.cpp"
    }

// We take the decision to reject or not to reject the null hypothesis H0
    for (i = 0; i <= (nblevel[0] - 1); i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
	//	if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }
    
// If applicable, we free the unused array of pointers

}

// We return
    return;
   
        
  }
  
}
