// Title: The Read-Cressie statistic for uniformity
// Ref. (book or article): Read, Timothy R. C. and Cressie, Noel A. C. (1988), Goodness-of-fit statistics for discrete multivariate data. 
//                         Springer Series in Statistics. Springer-Verlag, New York, 1988. xii+211 pp. ISBN: 0-387-96682-X (Reviewer: K. M. Lal Saxena)
// In fact, it is (5.5) page 459 of : Cressie, N. and Read, T.R.C., 1984, Multinomial goodness-of-fit tests. Journal of the Royal Statistic Society
// Series B, 46, 440–464. See also Marhuenda (2005) page 319


#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat72(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$2nI^{\\lambda}$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 1;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	paramstat[0] = 2.0 / 3.0;
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
    double lambda;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 1;
      lambda = 2.0 / 3.0;
      paramstat[0] = lambda;
    } else if (nbparamstat[0] == 1) {
      lambda = paramstat[0];
    } else {
      Rf_error("Number of parameters should be at most: 1");
    }
	
// If necessary, we check if some parameter values are out of parameter space

    
    if (n > 3) {
// Computation of the value of the test statistic
    void R_rsort (double* x, int n);
    double punif(double q, double min, double max, int lower_tail, int log_p);
    double *U;
    double *V;
    U = new double[n];
    V = new double[n + 1];
    double statI, sumI = 0.0;
	

    // generate vector U
    for (i = 0; i < n; i++) {
      U[i] = punif(x[i], 0.0, 1.0, 1, 0);
    }
    R_rsort(U, n); // We sort the data
	
    V[0] = U[0];
    V[n] = 1.0 - U[n - 1];	
    for (i = 1; i < n; i++) {
      V[i] = U[i] - U[i - 1];
    }
    
    // calculate statI
    for (i = 0; i < (n + 1); i++) {
      sumI = sumI + V[i] * (R_pow(((double)n + 1.0) * V[i], lambda) - 1.0);
    }
    
    statI = (double)(2 * n) * sumI / (lambda * (lambda + 1.0));
	
    statistic[0] = statI; // Here is the test statistic value
	

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue72.cpp"
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
    delete[] U;
    delete[] V;

}

// We return
    return;
   
        
  }
  
  
  
}
