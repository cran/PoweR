// Title: The HANSP statistic for the Laplace distribution
// Ref. (book or article): Tests of fit for the Laplace distribution based on correcting
// moments of entropy estimators. Hadi Alizadeh Noughabi and Sangun Park.
// JOURNAL OF STATISTICAL COMPUTATION AND SIMULATION, 2016
// VOL. 86, NO. 11, 2165–2181

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat84(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$AP_{l}$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 3;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	if ((int)paramstat[0] == 2) {
	  paramstat[0] = 2.0;
	  if ((int)paramstat[1] == 5) {
	    nom = "$AP_{z}^{(MLE)}$";
	    paramstat[1] = 5.0;
	  } else if ((int)paramstat[1] == 4) {
	    nom = "$AP_{a}^{(MLE)}$";
	    paramstat[1] = 4.0;
	  } else if((int)paramstat[1] == 3) {
	    nom = "$AP_{y}^{(MLE)}$";
	    paramstat[1] = 3.0;
	  } else if((int)paramstat[1] == 2) {
	    nom = "$AP_{e}^{(MLE)}$";
	    paramstat[1] = 2.0;
	  } else {
	    nom = "$AP_{v}^{(MLE)}$";
	    paramstat[1] = 1.0;
	  }
	} else {
	    paramstat[0] = 1.0;
	    if ((int)paramstat[1] == 5) {
	      nom = "$AP_{z}$";
	      paramstat[1] = 5.0;
	    } else if ((int)paramstat[1] == 4) {
	      nom = "$AP_{a}$";
	      paramstat[1] = 4.0;
	    } else if ((int)paramstat[1] == 3) {
	      nom = "$AP_{y}$";
	      paramstat[1] = 3.0;
	    } else if ((int)paramstat[1] == 2) {
	      nom = "$AP_{e}$";
	      paramstat[1] = 2.0;
	    } else {
	      nom = "$AP_{v}$";
	      paramstat[1] = 1.0;
	    }
	}
	if ((int)paramstat[0] == 1 && (int)paramstat[1] == 5) {
	  if (n == 0) paramstat[2] = 0.0;
	  if ((1 <= n) && (n <= 8)) paramstat[2] = 1.0;
	  if ((9 <= n) && (n <= 25)) paramstat[2] = 2.0;
	  if ((26 <= n) && (n <= 40)) paramstat[2] = 3.0;
	  if ((41 <= n) && (n <= 60)) paramstat[2] = 4.0;
	  if ((61 <= n) && (n <= 90)) paramstat[2] = 5.0;
	  if ((91 <= n) && (n <= 120)) paramstat[2] = 6.0;
	  if (121 <= n) paramstat[2] = round(R_pow(log((double)n), 1.15));
	} else {
	  if (n == 0) paramstat[2] = 0.0;
	  if ((1 <= n) && (n <= 8)) paramstat[2] = 1.0;
	  if ((9 <= n) && (n <= 15)) paramstat[2] = 2.0;
	  if ((16 <= n) && (n <= 25)) paramstat[2] = 4.0;
	  if ((26 <= n) && (n <= 40)) paramstat[2] = 5.0;
	  if ((41 <= n) && (n <= 60)) paramstat[2] = 6.0;
	  if ((61 <= n) && (n <= 90)) paramstat[2] = 7.0;
	  if ((91 <= n) && (n <= 120)) paramstat[2] = 8.0;
	  if (121 <= n) paramstat[2] = round(R_pow(log((double)n), 1.35));
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
    int group, version, m;
    // group is: First group or Second group of test statistics (1 or 2)
    // version is: v, e, y, a or z (1 to 5)
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 3;
      group = 1;
      version = 1;
      m = 1;
      paramstat[0] = 1.0;
      paramstat[1] = 1.0;
      paramstat[2] = 1.0;
    } else if (nbparamstat[0] == 1) {
      nbparamstat[0] = 3;
      group = (int)paramstat[0];
      version = 1;
      m = 1;
      paramstat[1] = 1.0;
      paramstat[2] = 1.0;
    } else if (nbparamstat[0] == 2) {
      nbparamstat[0] = 3;
      group = (int)paramstat[0];
      version = (int)paramstat[1];
      if (group == 1 && version == 5) {
	if ((1 <= n) && (n <= 8)) m = 1;
	if ((9 <= n) && (n <= 25)) m = 2;
	if ((26 <= n) && (n <= 40)) m = 3;
	if ((41 <= n) && (n <= 60)) m = 4;
	if ((61 <= n) && (n <= 90)) m = 5;
	if ((91 <= n) && (n <= 120)) m = 6;
	if (121 <= n) m = (int)round(R_pow(log((double)n), 1.15));
      } else {
	if ((1 <= n) && (n <= 8)) m = 1;
	if ((9 <= n) && (n <= 15)) m = 2;
	if ((16 <= n) && (n <= 25)) m = 4;
	if ((26 <= n) && (n <= 40)) m = 5;
	if ((41 <= n) && (n <= 60)) m = 6;
	if ((61 <= n) && (n <= 90)) m = 7;
	if ((91 <= n) && (n <= 120)) m = 8;
	if (121 <= n) m = (int)round(R_pow(log((double)n), 1.35));
      }
      paramstat[2] = (double)m;
    } else if (nbparamstat[0] == 3) {
      group = (int)paramstat[0];
      version = (int)paramstat[1];
      m = (int)paramstat[2];
      if (m == 0) {
	if (group == 1 && version == 5) {
	  if ((1 <= n) && (n <= 8)) m = 1;
	  if ((9 <= n) && (n <= 25)) m = 2;
	  if ((26 <= n) && (n <= 40)) m = 3;
	  if ((41 <= n) && (n <= 60)) m = 4;
	  if ((61 <= n) && (n <= 90)) m = 5;
	  if ((91 <= n) && (n <= 120)) m = 6;
	  if (121 <= n) m = (int)round(R_pow(log((double)n), 1.15));
	} else {
	  if ((1 <= n) && (n <= 8)) m = 1;
	  if ((9 <= n) && (n <= 15)) m = 2;
	  if ((16 <= n) && (n <= 25)) m = 4;
	  if ((26 <= n) && (n <= 40)) m = 5;
	  if ((41 <= n) && (n <= 60)) m = 6;
	  if ((61 <= n) && (n <= 90)) m = 7;
	  if ((91 <= n) && (n <= 120)) m = 8;
	  if (121 <= n) m = (int)round(R_pow(log((double)n), 1.35));
	}
      }
    } else {
      Rf_error("Number of parameters should be at most: 3");
    }

    // If necessary, we check if some parameter values are out of parameter space
    if ((group != 1) && (group != 2)) {
      Rf_warning("group should be 1 or 2 in stat84!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }
    if ((version != 1) && (version != 2) && (version != 3) && (version != 4) && (version != 5)) {
      Rf_warning("version should be 1, 2, 3, 4 or 5 in stat84!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }
    if (m > (n / 2)) {
      Rf_warning("m should be <= n/2 in stat84!\n");
      for (i = 0; i < n; i++) x[i] = R_NaN;
      return;
    }

    if (n > 3) {
// Computation of the value of the test statistic
      void R_rsort(double* x, int n);
      double temp1, temp2, temp;
      int k;
      double statval = 0.0;
      double Fhat84(int i, int n, double *x);
      double thetahatfunc84(double *beta, int n);
      
      R_rsort(x, n); 		// we sort the data
      //	    Rprintf(" %g ", Fhat84(21476, n, x)); Bug when n is too large due to interger overflow

      double bhat, muhat, tmp;
      if (group == 2) {
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
	tmp = 0.0;
	for (i = 0; i < n; i++) {
	  tmp = tmp + fabs(x[i] - muhat);
	}
	bhat = tmp / (double)n;
      }


      
      double *xi, *eta, *nu   = nullptr, *tau   = nullptr;
      if (group == 1) {
	
	xi = new double[n + 1];
	for (i = 1; i <= m; i++) {
	  temp = 0.0;
	  for (k = 2; k <= (i + m - 1); k++) temp = temp + x[k - 1];
	  xi[i - 1] = (temp + (double)(m + 2 - i) * x[0]) / (double)(2 * m);
	}
	for (i = (m + 1); i <= (n - m + 1); i++) {
	  temp = 0.0;
	  for (k = (i - m); k <= (i + m - 1); k++) temp = temp + x[k - 1];				
	  xi[i - 1] = temp / (double)(2 * m);
	}
	for (i = (n - m + 2); i <= (n + 1); i++) {
	  temp = 0.0;
	  for (k = (i - m); k <= (n - 1); k++) temp = temp + x[k - 1];
	  xi[i - 1] = (temp + (double)(i - n + m) * x[n - 1]) / (double)(2 * m);
	}
	
	if (version == 2) {
	  eta = new double[n + 1];
	  for (i = 1; i <= m; i++) {
	    temp = 0.0;
	    temp1 = x[0];
	    for (k = i; k <= m; k++) temp = temp + (x[m + k - 1] - temp1) / (double)(m + k - 1);
	    eta[i - 1] = xi[m] - temp;
	  }
	  for (i = (m + 1); i <= (n - m + 1); i++) eta[i - 1] = xi[i - 1];
	  for (i = (n - m + 2); i <= (n + 1); i++) {
	    temp = 0.0;
	    temp1 = x[n - 1];
	    for (k = (n - m + 2); k <= i; k++) temp = temp + (temp1 - x[k - m - 2]) / (double)(n + m - k + 1);
	    eta[i - 1] = xi[n - m] + temp;
	  }
	}

	if (version == 4) {
	  nu = new double[n + 1];
	  for (i = 1; i <= m; i++) {
	    temp = 0.0;
	    temp1 = x[0];
	    for (k = i; k <= m; k++) temp = temp + (x[m + k - 1] - temp1);
	    nu[i - 1] = xi[m] - temp / (double)m;
	  }
	  for (i = (m + 1); i <= (n - m + 1); i++) nu[i - 1] = xi[i - 1];
	  for (i = (n - m + 2); i <= (n + 1); i++) {
	    temp = 0.0;
	    temp1 = x[n - 1];
	    for (k = (n - m + 2); k <= i; k++) temp = temp + (temp1 - x[k - m - 2]);
	    nu[i - 1] = xi[n - m] + temp / (double)m;
	  }
	}

	if (version == 5) {
	  tau = new double[n + 1];
	  for (i = 1; i <= m; i++) {
	    temp = 0.0;
	    temp1 = x[0];
	    for (k = i; k <= m; k++) temp = temp + (x[m + k - 1] - temp1) / (double)k;
	    tau[i - 1] = xi[m] - temp;
	  }
	  for (i = (m + 1); i <= (n - m + 1); i++) tau[i - 1] = xi[i - 1];
	  for (i = (n - m + 2); i <= (n + 1); i++) {
	    temp = 0.0;
	    temp1 = x[n - 1];
	    for (k = (n - m + 2); k <= i; k++) temp = temp + (temp1 - x[k - m - 2]) / double(n - k + 2);
	    tau[i - 1] = xi[n - m] + temp;
	  }
	}
	
      } // End group == 1


      double Hv = 0.0;
      if (version == 1) {
	Hv = 0.0;
	for (i = 1; i <= n; i++) {
	  if ((i + m) <= n) temp1 = x[i + m - 1]; else temp1 = x[n - 1];
	  if ((i - m) >= 1) temp2 = x[i - m - 1]; else temp2 = x[0];
	  Hv = Hv + log((double)n * (temp1 - temp2) / ((double)(2 * m)));
	}
	Hv = Hv / (double)n;
      }
      
      double He, ci;
      if (version == 2) {
	He = 0.0;
	for (i = 1; i <= n; i++) {
	  if ((i + m) <= n) temp1 = x[i + m - 1]; else temp1 = x[n - 1];
	  if ((i - m) >= 1) temp2 = x[i - m - 1]; else temp2 = x[0];
	  if (i <= m) {
	    ci = 1.0 + (double)(i - 1) / (double)m;
	  } else {
	    if (i <= (n - m)) {
	      ci = 2.0;
	    } else {
	      ci = 1.0 + (double)(n - i) / (double)(m);
	    }
	  }
	  He = He + log((double)n * (temp1 - temp2) / (ci * (double)m));
	}
	He = He / (double)n;
      }    
	
      double Hy = 0.0, term;
      if (version == 3)  {
	Hy = 0.0;
	term = 0.0;
	for (i = 1; i <= n; i++) {
	  if ((i + m) <= n) {
	    temp1 = x[i + m - 1];
	    temp = Fhat84(i + m, n, x);
	  } else {
	    temp1 = x[n - 1];
	    temp = Fhat84(n, n, x);
	  }
	  if ((i - m) >= 1) {
	    temp2 = x[i - m - 1];
	    temp = temp - Fhat84(i - m, n, x);
	  } else {
	    temp2 = x[0];
	    temp = temp - Fhat84(1, n, x);
	  }
	  term = term + temp;
	  Hy = Hy + temp * log((temp1 - temp2) / temp);
	}
		
		Hy = Hy / term;

		
      }

      double *zeta  = nullptr;
      if (group == 1) {
      	if (version == 3) {
	  zeta = new double[n + 1];
	  for (i = 1; i <= (n + 1); i++) {
	    zeta[i - 1] = (double)(2 * m) * xi[i - 1] / term;
	  }
	}
      }

      double Ha = 0.0;
      int ai;
      if (version == 4) {
	Ha = 0.0;
	for (i = 1; i <= n; i++) {
	  if ((i + m) < n) {
	    temp1 = x[i + m - 1];
	  } else {
	    temp1 = x[n - 1];
	  }
	  if ((i - m) > 1) {
	    temp2 = x[i - m - 1];
	  } else {
	    temp2 = x[0];
	  }
	  if ((i <= m) || (i >= (n - m + 1))) ai = 1; else ai = 2;
	  Ha = Ha + log((double)n * (temp1 - temp2) / (double)(ai * m));
	}
	Ha = Ha / (double)n;
      }

      double Hz = 0.0, zi;
      if (version == 5)  {
	Hz = 0.0;
	for (i = 1; i <= n; i++) {
	  if ((i + m) < n) {
	    temp1 = x[i + m - 1];
	  } else {
	    temp1 = x[n - 1];
	  }
	  if ((i - m) > 1) {
	    temp2 = x[i - m - 1];
	  } else {
	    temp2 = x[0];
	  }
	  if (i <= m) {
	    zi = (double)i / (double)m;
	  } else if (i >= (n - m + 1)) {
	    zi = (double)(n - i + 1) / (double)m;
	  } else {
	    zi = 2.0;
	  }
	  Hz = Hz + log((double)n * (temp1 - temp2) / (zi * (double)m));
	}
	Hz = Hz / (double)n;
      }

      if (group == 1) {	
	if (version == 1) {
	  statval = log(2.0 * thetahatfunc84(xi, n)) + 1.0 - Hv; //APv
	} else if (version == 2) {
	  statval = log(2.0 * thetahatfunc84(eta, n)) + 1.0 - He; // APe
	} else if (version == 3) {
	  statval = log(2.0 * thetahatfunc84(zeta, n)) + 1.0 - Hy; // APy
	} else if (version == 4) {
	  statval = log(2.0 * thetahatfunc84(nu, n)) + 1.0 - Ha; // APa
	} else { // if (version == 5) {
	  statval = log(2.0 * thetahatfunc84(tau, n)) + 1.0 - Hz; // APz
	}
      } else if (group == 2) {
	if (version == 1) {
	  statval = log(2.0 * bhat) + 1.0 - Hv; // APvMLE
	} else if (version == 2) {
	  statval = log(2.0 * bhat) + 1.0 - He; // APeMLE
	} else if (version == 3) {
	  statval = log(2.0 * bhat) + 1.0 - Hy; // APyMLE
	} else if (version == 4) {
	  statval = log(2.0 * bhat) + 1.0 - Ha; // APaMLE
	} else { //if (version == 5) {
	  statval = log(2.0 * bhat) + 1.0 - Hz; // APzMLE
	}
      }

      statistic[0] = statval; // Here is the test statistic value
	

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue84.cpp"
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
    if (group == 1) {
      delete[] xi;
      if (version == 2) delete[] eta; 
      if (version == 3) delete[] zeta;
      if (version == 4) delete[] nu;
      if (version == 5) delete[] tau;
    }
}

// We return
    return;
   
        
  }

  double thetahatfunc84(double *beta, int n) {
    int i;
    double temp1, temp2, res;    
    if(n % 2 == 0) { // n even
      temp1 = 0.0;
      temp2 = 0.0;
      for (i = 1; i <= (n / 2); i++) temp1 = temp1 + 0.5 * (beta[i - 1] + beta[i]);
      for (i = ((n / 2) + 1); i <= n; i++) temp2 = temp2 + 0.5 * (beta[i - 1] + beta[i]);
      res = (temp2 - temp1) / (double)n;
    } else { // n odd
      temp1 = 0.0;
      temp2 = 0.0;
      for (i = 1; i <= ((n - 1) / 2); i++) temp1 = temp1 + 0.5 * (beta[i - 1] + beta[i]);
      for (i = (((n + 1)/ 2) + 1); i <= n; i++) temp2 = temp2 + 0.5 * (beta[i - 1] + beta[i]);
      res = (temp2 - temp1 + 0.25 * (beta[(n + 1) / 2] - beta[((n + 1) / 2 ) - 1])) / (double)n;
    }
    return res;
  }

  double Fhat84(int i, int n, double *x) {
    double rez = 0.0;
    double verysmall = R_pow(10.0, -50.0);
    if (i == 1) {
      rez = (double)n / (double)(n - 1);
      if (fabs(x[0] - x[1]) >= verysmall) rez = rez + (double)n / (double)(2 * n - 1);
    }
    if ((i >= 2) && (i <= (n - 1))) {
      rez = (double)(i * (n - 1) + 1) / (double)(n - 1);
      //      if (i == 21476) {Rprintf(" %g ", (double)(i * (n - 1) + 1)); Rprintf(" %g ", (double)(n - 1)); Rprintf(" (1): %g ", rez);}
      if (fabs(x[i] - x[i - 2]) >= verysmall) rez = rez + (x[i - 1] - x[i - 2]) / (x[i] - x[i - 2]);
    }
    if (i == n) {
      rez = (double)(n * n - 2 * n + 2) / (double)(n - 1);
      if (fabs(x[n - 3] - x[n - 1]) < verysmall) {
	//	rez = rez + 0.0;
      } else if ((fabs(x[n - 2] - x[n - 1]) < verysmall) && (fabs(x[n - 3] - x[n - 2]) > verysmall)) {
	rez = rez + 1.0;
      } else {
	rez = rez + 1.0 + (double)(n - 1) / (double)(2 * n - 1);
      }												  }
    return((double)(n - 1) * rez / (double)(n * (n + 1)));
  }

  // double muhatfunc(double *beta, int n) {
  //   double res;    
  //   if(n % 2 == 0) { // n even
  //     res = beta[n / 2];
  //   } else { // n odd
  //     res = 0.5 * (beta[(n + 1) / 2] + beta[(n + 1) / 2 - 1]);
  //   }
  //   return res;
  // }
  
}
