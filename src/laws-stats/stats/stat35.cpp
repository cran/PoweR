// Title: Statistique de test P_1 dans l'article 2 avec Alain et Alexandre ...
// Ref. (book or article): 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat35(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    if (alter[0] != 0 && alter[0] != 1 && alter[0] != 2) error("alter should be in {0,1,2}");

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$P_1$";
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
      for (i=j;i<50;i++) name[i][0] = space[0];
      return;
    }

    if (n>3) {
// Computation of the value of the test statistic
    double qnorm(double p, double mean, double sd, int lower_tail, int log_p);
    //    double pnorm(double x, double mu, double sigma, int lower_tail, int log_p);
    double *ychap;
    ychap = new double[n];
    double gamma=0.577215664901533, muchap=0.0, sigchap, sigchap2=0.0, deltachap=0.0, statvalue;
     

    for (i=0;i<=(n-1);i++) muchap = muchap + x[i];
    muchap = muchap/(double)n;

    for (i=1;i<=n;i++) sigchap2 = sigchap2 + R_pow(*(x+i-1)-muchap,2.0);
    sigchap2 = sigchap2/(double)n;
    sigchap = sqrt(sigchap2);
    for (i=0;i<=(n-1);i++) ychap[i] = (x[i]-muchap)/sigchap;
    for (i=0;i<=(n-1);i++) {
      if (ychap[i] != 0) {
	ychap[i] = R_pow(ychap[i],2.0) * log(R_pow(ychap[i],2.0));
      }
      else ychap[i] = 0.0;
    }
    for (i=0;i<=(n-1);i++) deltachap = deltachap + ychap[i];
    deltachap = deltachap/(double)n;

    statvalue = sqrt((double)n)*((deltachap-(2.0-log(2.0)-gamma))/sqrt(2.0*(3.0*R_pow(M_PI,2.0)/4.0-7.0)));

    statistic[0] = statvalue; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue35.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
      
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (alter[0] == 0) { if (statistic[0] > critvalR[i] || statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two-sided
	  } else if (alter[0] == 1) {  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0;  // less
	    } else { if (alter[0] == 2) {  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; } } // greater
    } else {
	  if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }
    
// If applicable, we free the unused array of pointers
    delete[] ychap;

}

// We return
    return;
    
        
  }
  
}
