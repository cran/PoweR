// Title: The Choi-Kim statistic for the Laplace distribution
// Ref. (book or article): Choi, B. and Kim, K. (2006), Testing goodness-of-fit for Laplace distribution based on maximum entropy,
//						   Statistics, Vol. 40, No. 6, pp. 517-531.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat52(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 4;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{m,n}^{E}$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	double mTE[] = {1.0, 2.0, 2.0,3.0, 3.0,4.0,4.0,5.0,2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
			2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
	if (n > 50) {
	  paramstat[0] = 2.0;	
	} else {
	  paramstat[0] = mTE[n - 4];
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
	// if parstat=NULL, par = 1
	// table 4 from choi2006 for $T_{m,n}^{E}$
    double mTE[] = {1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
		    2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    
    int m;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      if (n > 50) {
	m = 2;
	paramstat[0] = 2.0;	
      } else {
	m = (int)mTE[n - 4];
	paramstat[0] = mTE[n - 4];
      }
    } else if (nbparamstat[0] == 1) {
      m = (int)paramstat[0];
    } else {
     Rf_error("Number of parameters should be at most: 1");
    }

// If necessary, we check if some parameter values are out of parameter space
    if ((m > (n / 2)) || (m < 1)) {
      Rf_warning("m should be a positive integer smaller than n/2 in stat52!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }
    
    // Here m < (n/2)
    if (n > 3) {
// Computation of the value of the test statistic
      int k;
      void R_rsort(double* x, int n);
      double *Y;
      Y = new double[n];   	
      double statTEmn, tmp = 0.0, bhat, muhat = 0.0, prod, prod2;
      double temp1, temp2;
	
      // calculate mu^ and b^ by using the maximum likelihood estimators 
      // mu^ = the sample median
      // b^ = 1/n * \sum_{i=1}^{n} |xi - mu^|
      
      // calculate mu^
      R_rsort(x, n); 			// we sort the data from gensample
      if(n % 2 == 0) {		// check if n is divisible by 2
	muhat = (x[(n / 2) - 1] + x[n / 2]) / 2.0;
      } else {
	muhat = x[n / 2];
      }
	
      // calculate b^
      for (i = 0; i < n; i++) {
	Y[i] = x[i] - muhat;
	tmp = tmp + fabs(Y[i]);
      }
      bhat = tmp / (double)n;

      prod = 1.0;
      for (i = 1; i <= n; i++) {
	if ((i + m) > n) temp1 = Y[n - 1]; else temp1 = Y[i + m - 1];
	if ((i - m) < 1) temp2 = Y[0]; else temp2 = Y[i - m - 1];
	prod = prod * (temp1 - temp2);
      }
	
      prod2 = 1.0;
      for (k = m; k <= n; k++) {
	prod2 = prod2 * exp(1.0 / (double)k);
      }
	
      statTEmn = R_pow(prod, 1.0 / (double)(n - m)) * prod2 / bhat;
      
      statistic[0] = statTEmn; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue52.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i = 0; i <= (nblevel[0] - 1); i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value) Here we reject for small values of statistic!
      } else {
		if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
        }
    }
    
// If applicable, we free the unused array of pointers
    delete[] Y;

}

// We return
    return;
   
        
  }

  
}
