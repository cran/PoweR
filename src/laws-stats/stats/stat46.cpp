// Title: The Kuiper statistic for the Laplace distribution
// Ref. (book or article): Puig, P. and Stephens, M. A. (2000). Tests of fit for the Laplace distribution, with applications. Technometrics 42, 417-424.
// Now, one can also test for an APD(lambda) for any lambda.

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat46(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;

    int i, j = 0, n = xlen[0];
    if (getname[0] == 1) {    
      // Here, INDICATE the name of your statistic
      const char *nom = "$\\sqrt{n}V$";
      // Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 3;
      // Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	paramstat[0] = 0.5;
	paramstat[1] = 1.0;
	paramstat[2] = 1.0;		
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
    double theta1, theta2, lambda;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 3;
      paramstat[0] = 0.5;
      paramstat[1] = 1.0;
      paramstat[2] = 1.0;
      theta1 = 0.5;
      theta2 = 1.0;
      lambda = 1.0;
    } else if (nbparamstat[0] == 1) {
      nbparamstat[0] = 3;
      paramstat[1] = 1.0;
      paramstat[2] = 1.0;
      theta1 = paramstat[0];
      theta2 = 1.0;
      lambda = 1.0;
    } else if (nbparamstat[0] == 2) {
      nbparamstat[0] = 3;
      paramstat[2] = 1.0;
      theta1 = paramstat[0];
      theta2 = paramstat[1];
      lambda = 1.0;
    } else if (nbparamstat[0] == 3) {
      theta1 = paramstat[0];
      theta2 = paramstat[1];
      lambda = paramstat[2];
    } else {
      Rf_error("Number of parameters in stat46 should be at most: 3");
    }

  // If necessary, we check if some parameter values are out of parameter space
    if (lambda < 1.0 - 0.000000000000001) {
      Rf_warning("lambda should be >=1 in stat46!\n");
      for (i = 0; i < n; i++) statistic[0] = R_NaN;
      return;
    }
    if ((theta1 >= 1.0) || (theta1 <= 0)) {
      Rf_warning("theta1 should be in (0,1) in stat46!\n");
      for (i = 0; i < n; i++) statistic[0] = R_NaN;
      return;
    }
    if (theta2 <= 0) {
      Rf_warning("theta2 should be > 0 in stat46!\n");
      for (i = 0; i < n; i++) statistic[0] = R_NaN;
      return;
    }
    
    if (n > 3) {
// Computation of the value of the test statistic
      double pgamma(double q, double shape, double scale, int lowertail, int logp);     
      double myf46(double x, void *info);
      double R_zeroin2(			/* An estimate of the root */
		       double ax,				/* Left border | of the range	*/
		       double bx,				/* Right border| the root is seeked*/
		       double fa, double fb,		/* f(a), f(b) */
		       double (*f)(double x, void *info),	/* Function under investigation	*/
		       void *info,				/* Add'l info passed on to f	*/
		       double *Tol,			/* Acceptable tolerance		*/
		       int *Maxit);
    double *info;				/* Add'l info passed on to f	*/
    info = new double[n + 2];
    info[0] = lambda;
    info[1] = (double)n;
    for (i = 0; i < n; i++) info[i + 2] = x[i];
    double *Tol;			/* Acceptable tolerance		*/
    Tol = new double[1];
    Tol[0] = 0.000000000001;
    int *Maxit;
    Maxit = new int[1];
    Maxit[0] = 1000;
    double minx, maxx;
    minx = x[0];
    maxx = x[0];
    for (i = 0; i < n; i++) {
      if (x[i] < minx) minx = x[i];
      if (x[i] > maxx) maxx = x[i];
    }
    double fminx, fmaxx;
    fminx = myf46(minx, info);
    fmaxx = myf46(maxx, info);
      
    void R_rsort (double* x, int n);
    double plaplace(double y);
    double *Ka;
    Ka = new double[n];
    double statK, tmp = 0.0, tmp2, bhat, muhat, Dplus, Dminus;
      R_rsort(x, n); 		// we sort the data

	double deltatheta, y, tmp1;
	deltatheta = (2.0 * R_pow(theta1, theta2) * R_pow(1.0 - theta1, theta2)) / (R_pow(theta1, theta2) + R_pow(1.0 - theta1, theta2));

	if (fabs(lambda - 1.0) < 0.000000000000001) { // i.e. lambda = 1.0; Laplace case
    
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
	
      // generate vector Ka
      for (i = 0; i < n; i++) {
	      y = (x[i] - muhat) / bhat;
      if (y < 0.0) tmp1 = -y / theta1; else tmp1 = 0.0;
      if (y > 0.0) tmp2 = y / (1.0 - theta1); else tmp2 = 0.0;
      Ka[i] = theta1 * (1.0 - pgamma(deltatheta * R_pow(tmp1, theta2) / lambda, 1.0 / theta2, 1.0, 1, 0)) +
	(1.0 - theta1) * pgamma(deltatheta * R_pow(tmp2, theta2) / lambda, 1.0 / theta2, 1.0, 1, 0);

      //	Ka[i] = plaplace((x[i] - muhat) / bhat);
      }

      } else if (fabs(lambda - 2.0) < 0.000000000000001) {// i.e. lambda = 2.0; Normal case

	
	muhat = 0.0;
	for (i = 0; i < n; i++) {
	  muhat = muhat + x[i];
	}
	muhat = muhat / (double)n;
	for (i = 0; i < n; i++) {
	  tmp = tmp + R_pow(x[i] - muhat, 2.0);
	}
	bhat = sqrt(tmp / (double)n);
	
	for (i = 0; i < n; i++) {
	        y = (x[i] - muhat) / bhat;
      if (y < 0.0) tmp1 = -y / theta1; else tmp1 = 0.0;
      if (y > 0.0) tmp2 = y / (1.0 - theta1); else tmp2 = 0.0;
      Ka[i] = theta1 * (1.0 - pgamma(deltatheta * R_pow(tmp1, theta2) / lambda, 1.0 / theta2, 1.0, 1, 0)) +
	(1.0 - theta1) * pgamma(deltatheta * R_pow(tmp2, theta2) / lambda, 1.0 / theta2, 1.0, 1, 0);

      //	  Ka[i] = Rf_pnorm5((x[i] - muhat) / bhat, 0.0, 1.0, 1, 0);
	}

      } else {

	    double invlambda = 1.0 / lambda;

	// uniroot()
      muhat = R_zeroin2(minx,				/* Left border | of the range	*/
			maxx,				/* Right border| the root is seeked*/
			fminx, fmaxx,		/* f(a), f(b) */
			myf46,	/* Function under investigation	*/
			info,				/* Add'l info passed on to f	*/
			Tol,			/* Acceptable tolerance		*/
			Maxit);
    // calculate sighat^
    for (i = 0; i < n; i++) {
      tmp = tmp + R_pow(fabs(x[i] - muhat), lambda);
    }
    bhat = R_pow(tmp / (double)n, invlambda);

    for (i = 0; i < n; i++) {
      // See equ (2.3) in Desgagne et al. (2021)
      y = (x[i] - muhat) / bhat;
      if (y < 0.0) tmp1 = -y / theta1; else tmp1 = 0.0;
      if (y > 0.0) tmp2 = y / (1.0 - theta1); else tmp2 = 0.0;
      Ka[i] = theta1 * (1.0 - pgamma(deltatheta * R_pow(tmp1, theta2) / lambda, 1.0 / theta2, 1.0, 1, 0)) +
	(1.0 - theta1) * pgamma(deltatheta * R_pow(tmp2, theta2) / lambda, 1.0 / theta2, 1.0, 1, 0);
    }

      }
      
	
      // calculate statKS
      Dplus = 1.0 / (double)n - Ka[0];
      Dminus = Ka[0];
      for (i = 1; i < n; i++) {
	tmp2 = (double)(i + 1) / (double)n - Ka[i];
	if (tmp2 > Dplus) Dplus = tmp2;
	tmp2 = Ka[i] - ((double)i) / (double)n;
	if (tmp2 > Dminus) Dminus = tmp2;
      }
    	
      statK = sqrt((double)n) * (Dplus + Dminus);	
	
      statistic[0] = statK; // Here is the test statistic value
	

      if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
#include "pvalues/pvalue46.cpp"
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
      delete[] Ka;
 delete[] info;
 delete[] Tol;
 delete[] Maxit;
      
    }

    // We return
    return;
    
        
  }
  
      int sgn46(double v) {
    return (v > 0) - (v < 0);
  }

  double myf46(double x, void *info) {
    double lambda = (double)((double*)info)[0];
    int i, n = (int)((double*)info)[1];
    double tmp = 0.0;
    for (i = 0; i < n; i++) {
      tmp = tmp + R_pow(fabs((double)((double*)info)[i + 2] - x), lambda - 1.0) * (double)sgn46((double)((double*)info)[i + 2] - x);
    }
    return(tmp);
  }
  
  // In stat42.cpp, we already defined this function so no need to include it here.
  // The cumulative Laplace distribution function with \mu = 0 and \theta = 1
  // fabs = returns the absolute value of x (a negative value becomes positive, positive value is unchanged). 
  // double plaplace(double x) {
    // double temp = 0.5 * exp(-fabs(x));
    // return (x <= 0.0) ? temp : 1.0 - temp;
  // }
  
  
}
