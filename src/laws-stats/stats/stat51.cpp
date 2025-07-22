// Title: The Choi-Kim statistic for the Laplace distribution
// Ref. (book or article): Choi, B. and Kim, K. (2006), Testing goodness-of-fit for Laplace distribution based on maximum entropy,
//						   Statistics, Vol. 40, No. 6, pp. 517-531.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat51(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 4;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
      // Here, INDICATE the name of your statistic
      const char *nom = "$T_{m,n}^{VEC}$";  // was $T_{m,n}^{V}$ in JSS but now also computes $T_{m,n}^{E} or $T_{m,n}^{C}
      // Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 2;
      // Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	if ((int)paramstat[0] == 2) {
	  nom = "$T_{m,n}^{C}$";
	  paramstat[0] = 2.0;
	  double mTC[] = {1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 3.0, 2.0, 3.0, 3.0, 3.0,
			  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0,
			  4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0,
			  5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
	  if (n >= 51) {
	    paramstat[1] = ceil((double)n / 10.0);	
	  } else {
	    paramstat[1] = mTC[n - 4];
	  }
	} else if ((int)paramstat[0] == 3) {
	  nom = "$T_{m,n}^{E}$";
	  paramstat[0] = 3.0;
	  double mTE[] = {1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0};
	  if (n >= 12) {
	    paramstat[1] = 2.0;	
	  } else {
	    paramstat[1] = mTE[n - 4];
	  }
	} else {
	  nom = "$T_{m,n}^{V}$";
	  paramstat[0] = 1.0;
	  double mTV[] = {1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
			  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0,
			  4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                  5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 6.0, 6.0, 6.0, 6.0};
	  if (n >= 51) {
	    paramstat[1] = ceil((double)n / 10.0);	
	  } else {
	    paramstat[1] = mTV[n - 4];
	  }
	}
	if (n == 0) paramstat[1] = 0.0;
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
    // Table 4 from choi2006 for $T_{m,n}^{V}$
    double mTV[] = {1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
		    3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0,
		    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                  5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 6.0, 6.0, 6.0, 6.0};
    double mTC[] = {1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 3.0, 2.0, 3.0, 3.0, 3.0,
		    3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0,
		    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0,
                  5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
    double mTE[] = {1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0};
    
    int version, m = 0;
    // version = 1 for CKv; 2 for CKc; and 3 for CKe
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 2;
      version = 1;
      paramstat[0] = 1.0;
      if (n >= 51) {
	m = (int)ceil((double)n / 10.0);
	paramstat[1] = ceil((double)n / 10.0);	
      } else {
	m = (int)mTV[n - 4];
        paramstat[1] = mTV[n - 4];
      }
    } else if (nbparamstat[0] == 1) {
      nbparamstat[0] = 2;
      version = (int)paramstat[0];
      if (version == 1) {
	if (n >= 51) {
	  m = (int)ceil((double)n / 10.0);
	  paramstat[1] = ceil((double)n / 10.0);	
	} else {
	  m = (int)mTV[n - 4];
	  paramstat[1] = mTV[n - 4];
	}
      } else if (version == 2) {
	if (n >= 51) {
	  m = (int)ceil((double)n / 10.0);
	  paramstat[1] = ceil((double)n / 10.0);	
	} else {
	  m = (int)mTC[n - 4];
	  paramstat[1] = mTC[n - 4];
	}
      } else if (version == 3) {
	if (n >= 12) {
	  m = 2;
	  paramstat[1] = 2.0;	
	} else {
	  m = (int)mTE[n - 4];
	  paramstat[1] = mTE[n - 4];
	}
      }
    } else if (nbparamstat[0] == 2) {
      version = (int)paramstat[0];
      m = (int)paramstat[1];
      if (m == 0) {
	if (version == 1) {
	  if (n >= 51) {
	    m = (int)ceil((double)n / 10.0);
	  } else {
	    m = (int)mTV[n - 4];
	  }
	} else if (version == 2) {
	  if (n >= 51) {
	    m = (int)ceil((double)n / 10.0);
	  } else {
	    m = (int)mTC[n - 4];
	  }
	} else if (version == 3) {
	  if (n >= 12) {
	    m = 2;
	  } else {
	    m = (int)mTE[n - 4];
	  }
	}
      }
    } else {
      Rf_error("Number of parameters should be at most: 2");
    }

// If necessary, we check if some parameter values are out of parameter space
    if ((version != 1) && (version != 2) && (version != 3)) {
      Rf_warning("version should be 1, 2 or 3 in stat51!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }
    if ((m > (n / 2)) || (m < 1)) {
      Rf_warning("m should be a positive integer smaller than n/2 in stat51!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }
	

    // Here m < (n/2)
    if (n > 3 && m < n / 2) {
      // Computation of the value of the test statistic
      void R_rsort(double* x, int n);
      int k, idx;
      double tmp = 0.0, bhat, muhat, prod = 1.0;
      double temp1, temp2, statval;
      
      // calculate mu^ and b^ by using the maximum likelihood estimators 
      // mu^ = the sample median
      // b^ = 1/n * \sum_{i=1}^{n} |xi - mu^|
      
      // calculate mu^
      R_rsort(x, n); 		// we sort the data
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

      if (version == 1) {
	prod = 0.0;
	for (i = 1; i <= n; i++) {
	  if ((i + m) >= n) temp1 = x[n - 1]; else temp1 = x[i + m - 1];
	  if ((i - m) <= 1) temp2 = x[0]; else temp2 = x[i - m - 1];
	  prod = prod + log(temp1 - temp2);
	}
	statval = (double)n * exp(prod / (double)n) / ((double)(2 * m) * bhat); // CKv
      }

      if (version == 2) {
	prod = 0.0;
	for (i = 1; i <= n; i++) {
	  tmp = 0.0;
	  for (k = (i - m); k <= (i + m); k++) {
	    if (k < 1) idx = 0; else if (k >= n) idx = n - 1; else idx = k - 1;
	    tmp = tmp + x[idx];
	  }
	  tmp = tmp / (double)(2 * m + 1);
	  temp1 = 0.0; temp2 = 0.0;
	  for (j = (i - m); j <= (i + m); j++) {
	    if (j < 1) idx = 0; else if (j >= n) idx = n - 1; else idx = j - 1;		    
	    temp1 = temp1 + (double)(j - i) * (x[idx] - tmp) / (double)n;
	    temp2 = temp2 + R_pow(x[idx] - tmp, 2.0);
	  }
	  prod = prod + log(temp1 / temp2);
	}
	statval = exp(-prod / (double)n) / bhat; // CKc
      }
      
      if (version == 3) {
	prod = 0.0;
	for (i = 1; i <= (n - m); i++) {
	  if ((i + m) >= n) temp1 = x[n - 1]; else temp1 = x[i + m - 1];
	  prod = prod + log(temp1 - x[i - 1]);
	}
	statval = exp(prod / (double)(n - m));
	tmp = 0.0;
	for (k = m; k <= n; k++) tmp = tmp + 1.0 / (double)k;
	tmp = exp(tmp) / bhat;
	statval = statval * tmp; // CKe
      }



      
      statistic[0] = statval; // Here is the test statistic value
      
      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue51.cpp"
      }
      
      // We take the decision to reject or not to reject the null hypothesis H0
      for (i = 0; i < nblevel[0]; i++) {
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value) Here we reject for small values of statistic!
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


