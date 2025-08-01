// Title: The Morales statistic for uniformity
// Ref. (book or article): Morales, D., Pardo, L., Pardo, M. C. and Vajda, I. (2003), Limit laws for disparities of spacings,
//						   Journal of Nonparametric Statistics, 15(3), 325-342.
// Y. Marhuenda, D. Morales & M. C. Pardo, (2005), A comparison of uniformity tests, 
// A Journal of Theoretical and Applied Statistics, 39:4, 315-327


#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat78(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$D_{n,m}(\\phi_\\lambda)$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 2;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
      paramstat[0] = 0.0;
      paramstat[1] = 2.0;
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
    int m;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 2;
      lambda = 0.0;
      m = 2;
      paramstat[0] = 0.0;
      paramstat[1] = 2.0;
    } else if (nbparamstat[0] == 1) {
      nbparamstat[0] = 2;
      lambda = paramstat[0];
      m = 2;
      paramstat[1] = 2.0;
    } else if (nbparamstat[0] == 2) {
      lambda = paramstat[0];
      m = (int)paramstat[1];
    } else {
      Rf_error("Number of parameters should be at most: 2");
    }

    // If necessary, we check if some parameter values are out of parameter space
    if (m > (n + 1)) {
      Rf_warning("m should be <= n+1 in stat78!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

    if (n>3) {
// Computation of the value of the test statistic
    void R_rsort (double* x, int n);
    double punif(double q, double min, double max, int lower_tail, int log_p);
    double *U;
    double *Gnm;
    U = new double[n];
    Gnm = new double[n];
    double temp, statDn;
	

    // generate vector U
    for (i = 0; i < n; i++) {
      U[i] = punif(x[i], 0.0, 1.0, 1, 0);
    }
    R_rsort(U, n); // We sort the data
	
    // generate vector Gnm
    for (i = 0; i < n; i++) {
      if (i < (n - m)) {
	Gnm[i] = U[i + m] - U[i];
      } else {
	Gnm[i] = 1.0 + U[i + m - n] - U[i];
      }
    }
	
    // calculate statDn
    if (fabs(lambda) < 0.0000000000000001) { // lambda = 0.0
      statDn = 0.0;
      for (i = 0; i < n; i++) {
	temp = (double)n * Gnm[i] / (double)m;
	if (Gnm[i] > 0.0000000000000001) statDn = statDn + temp * log(temp);
      }
    } else if (fabs(lambda + 1.0) < 0.0000000000000001) { // lambda = -1.0
      statDn = 0.0;
      for (i = 0; i < n; i++) {
	statDn = statDn - log((double)n * Gnm[i] / (double)m);
      }
    } else {
      statDn = 0.0;
      for (i = 0; i < n; i++) {
	statDn = statDn + (R_pow((double)n * Gnm[i] / (double)m, lambda + 1.0) - 1.0) / (lambda * (lambda + 1.0));
      }
    }
    
    statistic[0] = statDn; // Here is the test statistic value
	

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue78.cpp"
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
    delete[] U;
    delete[] Gnm;

}

// We return
    return;
   
        
  }
  
  
  
}
