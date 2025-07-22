// Title: Subramanian / Dixit test
// Ref. (book or article): L. SUBRAMANIAN AND V. U. DIXIT, Tests for the skewness parameter of two-piece double
// exponential distribution in the presence of nuisance parameters, STATISTICS, 2017

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat108(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 0;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {
      // Here, INDICATE the name of your statistic
      const char *nom = "$SD$";
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
      double hsm108(double *x, int nx);
      void R_rsort (double* x, int n);

      int n1 = 0;
      double u, v, theta;
      
      R_rsort(x, n); 		// we sort the data

      theta = hsm108(x, n);
      
      for (i = 0; i < n; i++) {
	if (x[i] <= theta) n1 = n1 + 1;
      }
     
      u = (double)n1 * x[n1 - 1];
      for (i = 0; i < n1; i++) {
	u = u - x[i];
      }

      if (n1 < n) {
	v = -(double)(n - n1) * x[n1];
	for (i = n1; i < n; i++) {
	  v = v + x[i];
	}
      } else {
	v = 0.0;
      }

      if (u > 0 || v > 0) {
	statistic[0] = u / (u + v); // Here is the test statistic value
      } else {
	statistic[0] = 0.5;
      }
      
      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue108.cpp"
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

  double hsm108(double *x, int nx) {
    // x is already sorted
    int j;
    double *xtmp;
    xtmp = new double[nx];
    for (j = 0; j < nx; j++) {
      xtmp[j] = x[j];
    }
    int k = nx / 2 + (nx % 2 != 0) - 1; // ceiling(nx / 2.0) - 1;
    int nxmk;
    nxmk = nx - k;
    double *inf, *sup, *diffs;
    inf = new double[nxmk];
    sup = new double[nxmk];
    diffs = new double[nxmk];

    double min = xtmp[nx - 1] - xtmp[0];

    int i;

    double mean;
    int sum, nsum;

    double z, M;
    
    while (nx >= 4) {

      min =  xtmp[nx - 1] - xtmp[0];
      for (j = 0; j < nxmk; j++) {
	inf[j] = xtmp[j];
	sup[j] = xtmp[k + j];
	diffs[j] = sup[j] - inf[j];
	if (diffs[j] < min) min = diffs[j];
      }

      if (nxmk > 1) {
	sum = 0; nsum = 0;
	for (j = 0; j < nxmk; j++) {
	  if (fabs(diffs[j] - min) < 0.000000000000001) {
	    sum = sum + j + 1;
	    nsum = nsum + 1;
	  }
	}
	i = (int)((double)sum / (double)nsum);
      } else {
	i = 1;
      }

      if (fabs(diffs[i - 1]) < 0.000000000000001) {
	nx = 1;
	xtmp[0] = xtmp[i - 1];
      } else {
	nx = k + 1;
	for (j = 0; j < nx; j++) {
	  xtmp[j] = xtmp[i + j - 1];
	}
      }
      k  = nx / 2 + (nx % 2 != 0) - 1; // ceiling(nx / 2.0) - 1;
      nxmk = nx - k;
    }
      
    if (nx == 3) {
      z = 2.0 * xtmp[1] - xtmp[0] - xtmp[2];
      if (z < 0) {
	M = 0.5 * (xtmp[0] + xtmp[1]);
      } else if (z > 0) {
	M = 0.5 * (xtmp[1] + xtmp[2]);
      } else {
	M = xtmp[1];
      }
    } else {
      mean = 0.0;
      for (j = 0; j < nx; j++) {
	mean = mean + xtmp[j];
      }
      M = mean / (double)nx;
    }
      
    delete[] inf;
    delete[] sup;
    delete[] diffs;
    delete[] xtmp;
    
    return(M); 
  }
  
}
