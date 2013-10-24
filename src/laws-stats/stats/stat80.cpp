// Title: The Marhuenda statistic for uniformity
// Ref. (book or article): M.A. Marhuenda, Y. Marhuenda, D. Morales, (2005), "Uniformity tests under quantile categorization", 
//						   Kybernetes, Vol. 34 Iss: 6, pp.888 - 901.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat80(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{n,m}^{\\lambda}$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 2;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
      paramstat[0] = 1.0;
      paramstat[1] = 5.0;
     }
// The following 7 lines should NOT be modified
      const char *space = " ";
      while (nom[j] != '\0') {
	name[j][0] = nom[j];
	j++;
      }
      for (i=j;i<50;i++) name[i][0] = space[0];
      return;
    }

// Initialization of the parameters
    double lambda, m;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 2;
      lambda = 1.0;
	  m = 5.0;
      paramstat[0] = 1.0;
      paramstat[1] = 5.0;
    } else if (nbparamstat[0] == 1) {
      nbparamstat[0] = 2;
      lambda = paramstat[0];
      m = 5.0;
      paramstat[1] = 5.0;
    } else if (nbparamstat[0] == 2) {
      lambda = paramstat[0];
      m = paramstat[1];
    } else {
      return;
    }
	

    if (n>3) {
// Computation of the value of the test statistic
    void R_rsort (double* x, int n);
    double punif(double q, double min, double max, int lower_tail, int log_p);
    double *U;
	double *Z1;
    double *Z2;
	U = new double[n];
	Z1 = new double[n];
	Z2 = new double[n];
	double statTn, sumTn=0.0, pni;
	

	// generate vector U
	for (i=0;i<n;i++) {
	  U[i] = punif(x[i],0.0,1.0,1,0);
	}
    R_rsort(U,n); // We sort the data
	
	    
    for (i=1; i<=((int)m); i++) {
      pni = U[(i*n)/(int)m - 1] - U[((i-1)*n)/(int)m - 1];
      sumTn = sumTn + pni*(R_pow(m*pni,lambda) - 1.0);
    }
	

	statTn = 2.0*(double)n*sumTn/(lambda*(lambda+1.0));
	
	
    statistic[0] = statTn; // Here is the test statistic value
	

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue80.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
      } else {
		  if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
        }
    }
    
// If applicable, we free the unused array of pointers
    delete[] U;
	delete[] Z1;
	delete[] Z2;

}

// We return
    return;
   
        
  }
  
  
  
}
