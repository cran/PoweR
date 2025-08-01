// Title: The AJ statistic for the Laplace distribution
// Ref. (book or article): Informational energy-based goodness-of-fit test for
// Laplace distribution
// Hadi Alizadeh Noughabi and Jalil Jarrahiferiz
// Int. J. Information and Decision Sciences, Vol. 11, No. 3, (2019)


#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat102(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$AJ$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	if (n == 0) paramstat[0] = 0.0;
	if ((1 <= n) && (n <= 8)) paramstat[0] = 2.0;
	if ((9 <= n) && (n <= 15)) paramstat[0] = 3.0;
	if ((16 <= n) && (n <= 25)) paramstat[0] = 5.0;
	if ((26 <= n) && (n <= 35)) paramstat[0] = 6.0;
	if ((36 <= n) && (n <= 45)) paramstat[0] = 7.0;
	if ((46 <= n) && (n <= 60)) paramstat[0] = 8.0;
	if ((61 <= n) && (n <= 90)) paramstat[0] = 9.0;
	if ((91 <= n) && (n <= 120)) paramstat[0] = 10.0;
	if (n >= 121) paramstat[0] = round(R_pow(log((double)n), 1.5));
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
    int m;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      if ((1 <= n) && (n <= 8)) m = 2;
      if ((9 <= n) && (n <= 15)) m = 3;
      if ((16 <= n) && (n <= 25)) m = 5;
      if ((26 <= n) && (n <= 35)) m = 6;
      if ((36 <= n) && (n <= 45)) m = 7;
      if ((46 <= n) && (n <= 60)) m = 8;
      if ((61 <= n) && (n <= 90)) m = 9;
      if ((91 <= n) && (n <= 120)) m = 10;
      if (n >= 121) m = (int)round(R_pow(log((double)n), 1.5));
      paramstat[0] = (double)m;
    } else if (nbparamstat[0] == 1) {
      m = (int)paramstat[0];
      if (m == 0) {
	if ((1 <= n) && (n <= 8)) m = 2;
	if ((9 <= n) && (n <= 15)) m = 3;
	if ((16 <= n) && (n <= 25)) m = 5;
	if ((26 <= n) && (n <= 35)) m = 6;
	if ((36 <= n) && (n <= 45)) m = 7;
	if ((46 <= n) && (n <= 60)) m = 8;
	if ((61 <= n) && (n <= 90)) m = 9;
	if ((91 <= n) && (n <= 120)) m = 10;
	if (n >= 121) m = (int)round(R_pow(log((double)n), 1.5));
      }
    } else {
      Rf_error("Number of parameters should be at most: 1");
    }

// If necessary, we check if some parameter values are out of parameter space
    if (m > (n / 2)) {
      Rf_warning("m should be <= n/2 in stat102!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

    
    if (n > 3) {
// Computation of the value of the test statistic
      void R_rsort(double* x, int n);
      double plaplace(double y);
      double temp1, temp2;
      
      R_rsort(x, n); 		// we sort the data

      double muhat, bhat, tmp = 0.0;
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
      for (i = 0; i < n; i++) {
	tmp = tmp + fabs(x[i] - muhat);
      }
      bhat = tmp / (double)n;


      tmp = 0.0;
      for (j = 1; j <= n; j++) {
	if ((j + m) <= n) temp1 = x[j + m - 1]; else temp1 = x[n - 1];
	if ((j - m) >= 1) temp2 = x[j - m - 1]; else temp2 = x[0];
	temp1 = plaplace((temp1 - muhat) / bhat);
	temp2 = plaplace((temp2 - muhat) / bhat);
	tmp = tmp + 1.0 / (temp1 - temp2);
      }
      tmp = ((double)(2 * m)) * tmp / (double)(n * n);
      
      statistic[0] = tmp; // Here is the test statistic value
      

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue102.cpp"
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

}

// We return
    return;
   
        
  }
 
}
