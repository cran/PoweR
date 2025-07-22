// Title: The Desgagne-Lafaye de Micheaux-Ouimet statistic for the Laplace distribution
// Ref. (book or article):

#include <R.h>
#include "Rmath.h"
#include <assert.h>
#define EPSILON DBL_EPSILON

extern "C" {

  void stat106(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    if (alter[0] != 0 && alter[0] != 1 && alter[0] != 2 && alter[0] != 3) Rf_error("alter should be in {0,1,2,3}");

    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {    
// Here, INDICATE the name of your statistic
      const char *nom = "$DLO(\\lambda,version)$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 2;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).
	paramstat[0] = 1.0;
	paramstat[1] = 1.0;
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
    int version;
    if (nbparamstat[0] == 0) {
      nbparamstat[0] = 2;
      paramstat[0] = 1.0;
      paramstat[1] = 1.0;
      lambda = 1.0;
      version = 1;
    } else if (nbparamstat[0] == 1) {
      nbparamstat[0] = 2;
      lambda = paramstat[0];
      paramstat[1] = 1.0;
      version = 1;
    } else if (nbparamstat[0] == 2) {
      lambda = paramstat[0];
      version = (int)paramstat[1];
    } else {
      Rf_error("Number of parameters in stat106 should be at most: 2");
    }

  // If necessary, we check if some parameter values are out of parameter space
    if (lambda < 1.0 - 0.000000000000001) {
      Rf_warning("lambda should be >=1 in stat106!\n");
      for (i = 0; i < n; i++) statistic[0] = R_NaN;
      return;
    }
    if ((version != 1) && (version != 2) && (version != 3) && (version != 4) && (version != 5) && (version != 6)) {
      Rf_warning("version should be in {1,2,3,4,5,6} in stat106!\n");
      for (i = 0; i < n; i++) statistic[0] = R_NaN;
      return;
    }

    
    if (n>3) {
// Computation of the value of the test statistic
      double myf106(double x, void *info);
      double R_zeroin2(			/* An estimate of the root */
		       double ax,				/* Left border | of the range	*/
		       double bx,				/* Right border| the root is seeked*/
		       double fa, double fb,		/* f(a), f(b) */
		       double (*f)(double x, void *info),	/* Function under investigation	*/
		       void *info,				/* Add'l info passed on to f	*/
		       double *Tol,			/* Acceptable tolerance		*/
		       int *Maxit);
      
    void R_rsort(double* x, int n);
    //    double Rf_pchisq(double q, double df, int lower_tail, int log_p)
      //    double Rf_pnorm5(double q, double mu, double sigma, int lower_tail,int log_p);
    double Rf_gammafn (double x);	// gamma function in R extensions
    double Rf_digamma(double x);	// digamma function in R extensions
    double Rf_trigamma (double x);	// trigamma function in R extensions
    int sgn106(double v);

    double invlambda = 1.0 / lambda;
    double beta = 1.0 + invlambda;
    double gamma = -Rf_digamma(1.0);
    double onemgamma = 1.0 - gamma;
    double tenmfifteen = 0.000000000000001;
    bool islambdaone = fabs(lambda - 1.0) < tenmfifteen;
    bool islambdatwo = fabs(lambda - 2.0) < tenmfifteen;
    bool islambdaonephalf = abs(lambda - 1.5) < tenmfifteen;
    bool islambdatwophalf = fabs(lambda - 2.5) < tenmfifteen;
    bool islambdathree = abs(lambda - 3.0) < tenmfifteen;
    double onesixteenth = 1.0 / 16.0;
    double pisquare = R_pow(M_PI, 2.0);
    double ctepi = pisquare / 3.0 - 3.0;
    double ctepibis = 3.0 - 8.0 / M_PI;
    double ctepibisbis = (3.0 * pisquare - 28.0) / 8.0;
    double ctegamma = 0.5 * (2.0 - log(2.0) - gamma);
    double ctebeta = invlambda * beta * Rf_trigamma(beta) - invlambda;
    double ctelambda = invlambda * (lambda + log(lambda) + Rf_digamma(invlambda));
    double ctelambdabis = 1.0 + lambda - R_pow(lambda, 2.0) / (Rf_gammafn(2.0 - invlambda) * Rf_gammafn(invlambda));
    double sqrtn = sqrt((double)n);
    
    double mean, muhat;
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
    fminx = myf106(minx, info);
    fmaxx = myf106(maxx, info);
    
      R_rsort(x, n); 			// we sort the data from gensample
    if (islambdaone) { // lambda = 1
      // calculate mu^
      if(n % 2 == 0) {		// check if n is divisible by 2
	muhat = (x[(n / 2) - 1] + x[n / 2]) / 2.0;
      } else {
	muhat = x[((n + 1) / 2) - 1];
      }
    } else if (islambdatwo) { // lambda = 2
      mean = 0.0;
      for (i = 0; i < n; i++) {
	mean = mean + x[i];
      }
      muhat = mean / (double)n;
  } else {
	// uniroot()
      muhat = R_zeroin2(minx,				/* Left border | of the range	*/
			maxx,				/* Right border| the root is seeked*/
			fminx, fmaxx,		/* f(a), f(b) */
			myf106,	/* Function under investigation	*/
			info,				/* Add'l info passed on to f	*/
			Tol,			/* Acceptable tolerance		*/
			Maxit);
  }
    // calculate sighat^
    double tmp = 0.0, sighat;
    for (i = 0; i < n; i++) {
      tmp = tmp + R_pow(fabs(x[i] - muhat), lambda);
    }
    sighat = R_pow(tmp / (double)n, invlambda);


    double Blambda = 0.0, Klambda = 0.0, y, absy, plafond;
    // Equ. (4.14)
    plafond = R_pow(10.0, 300.0 / lambda);
    for (i = 0; i < n; i++) {
      y = (x[i] - muhat) / sighat;
      absy = fabs(y);
      if (absy > plafond) absy = plafond;
      tmp = R_pow(absy, lambda);
      Blambda = Blambda + tmp * sgn106(y);
      if (absy > 0.000000000000001) Klambda = Klambda + tmp * log(absy);
    }
    Blambda = Blambda / (double)n;
    Klambda = Klambda / (double)n;

    double denom1x, num2x, denom2x, term1x, term2x, XlambdaAPDstar;
    if (islambdaone) { // lambda = 1
      // Equ. (4.16)
      denom1x = 1.0;
      num2x = onemgamma;
      denom2x = ctepi;
    } else if (islambdatwo) { // lambda = 2
      // Equ. (4.17)
      denom1x = ctepibis;
      num2x = 0.5 * (2.0 - log(2.0) - gamma);
      denom2x = ctepibisbis;
    } else {
      denom1x = ctelambdabis;
      num2x = ctelambda;
      denom2x = invlambda * (beta * Rf_trigamma(beta) - 1.0);
  }
    // Equ. (4.15)
    term1x = (double)n * R_pow(Blambda, 2.0) / denom1x;
    term2x = (double)n * R_pow(Klambda - num2x, 2.0) / denom2x;
    XlambdaAPDstar = term1x + term2x; // X_{\lambda}^{APD*}


    double c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0, alpha1 = 0.0, alpha2 = 0.0, alpha3 = 0.0, alpha4 = 0.0;
    if (islambdaone) { // lambda = 1
      if (n % 2 == 0) { // n even
	alpha1 = 1.06;
	c1 = -1.856;
	alpha2 = 1.01;
	c2 = -0.422;
	alpha3 = 0.92;
	c3 = -1.950;
	alpha4 = 2.3;
	c4 = 39.349;
      } else { // n odd
	alpha1 = 1.03;
	c1 = -0.281;
	alpha2 = 0.86;
	c2 = -0.198;
	alpha3 = 1.04;
	c3 = -3.827;
	alpha4 = 1.0;
	c4 = 0;
      }
    } else if (islambdatwo) { // lambda = 2
      alpha1 = 0.99;
      c1 = -1.890;
      alpha2 = 1.00;
      c2 = -0.788;
      alpha3 = 1.05;
      c3 = -9.327;
      alpha4 = 1.4;
      c4 = 14.208;
    } else if (islambdaonephalf) { // lambda = 1.5
      alpha1 = 0.99;
      c1 = -0.952;
      alpha2 = 0.99;
      c2 = -0.637;
      alpha3 = 0.55;
      c3 = -3.488;
      alpha4 = 0.5;
      c4 = 2.434;
    } else if (islambdatwophalf) { // lambda = 2.5
      alpha1 = 0.99;
      c1 = -2.981;
      alpha2 = 0.99;
      c2 = -0.844;
      alpha3 = 1.10;
      c3 = -23.104;
      alpha4 = 1.3;
      c4 = 30.028;
    } else if (islambdathree) { // lambda = 3
      alpha1 = 0.97;
      c1 = -3.855;
      alpha2 = 0.98;
      c2 = -0.880;
      alpha3 = 1.14;
      c3 = -95.743;
      alpha4 = 1.2;
      c4 = 103.871;
    }


    double NKlambda, VarBlambda, VarNKlambdaonefourth, ENKlambdaonefourth, tmplambda, ZBlambda, ZNKlambdaonefourth, XlambdaAPD;
    // Equ. (5.1)
    NKlambda = Klambda - 0.5 * lambda * R_pow(Blambda, 2.0);
    if (NKlambda < 0.0) NKlambda = 0.0;
    if (islambdaone) { // lambda = 1
      if (n % 2 == 0) { // n even
	// Equ. (5.9) and (5.10)
	VarBlambda = (1.0 - 1.856 / R_pow((double)n, 1.06)) / (double)n;
	VarNKlambdaonefourth = (onesixteenth * R_pow(onemgamma, -1.5) * ctepi * (1.0 - 1.950 / R_pow((double)n, 0.92) + 39.349 / R_pow((double)n, 2.3))) / (double)n;
	ENKlambdaonefourth = R_pow(onemgamma, 0.25) * (1.0 - 0.422 / R_pow((double)n, 1.01));
      } else { // n odd
	// Equ. (5.11) and (5.12)
	VarBlambda = (1.0 - 0.281 / R_pow((double)n, 1.03)) / (double)n;
	VarNKlambdaonefourth = (onesixteenth * R_pow(onemgamma, -1.5) * ctepi * (1.0 - 3.827 / R_pow((double)n, 1.04))) / (double)n;
	ENKlambdaonefourth = R_pow(onemgamma, 0.25) * (1.0 - 0.198 / R_pow((double)n, 0.86));
      }
    } else if (islambdatwo) { // lambda = 2
      // Equ. (5.13) and (5.14)
      VarBlambda = ctepibis * (1.0 - 1.890 / R_pow((double)n, 0.99)) / (double)n;
      VarNKlambdaonefourth = (onesixteenth * R_pow(ctegamma, -1.5) * ctepibisbis * (1.0 - 9.327 / R_pow((double)n, 1.05) + 14.208 / R_pow((double)n, 1.4))) / (double)n;
      ENKlambdaonefourth = R_pow(ctegamma, 0.25) * (1.0 - 0.788 / (double)n);
    } else {
      // Equ. (5.4) and (5.5) and (5.6)
      VarBlambda = ctelambdabis * (1.0 + c1 / R_pow((double)n, alpha1)) / (double)n;
      tmplambda = ctelambda;
      ENKlambdaonefourth = R_pow(tmplambda, 0.25) * (1.0 + c2 / R_pow((double)n, alpha2));
      VarNKlambdaonefourth = onesixteenth * R_pow(tmplambda, -1.5) * invlambda * (beta * Rf_trigamma(beta) - 1.0) * (1.0 + c3 / R_pow((double)n, alpha3) + c4 / R_pow((double)n, alpha4)) / (double)n;
    }



    if (islambdaone || islambdatwo || islambdaonephalf || islambdatwophalf || islambdathree) {
      // Equ. (5.7) and (5.8) and (5.15)
      ZBlambda = Blambda / sqrt(VarBlambda);
      ZNKlambdaonefourth = (R_pow(NKlambda, 0.25) - ENKlambdaonefourth) / sqrt(VarNKlambdaonefourth);
      XlambdaAPD = R_pow(ZBlambda, 2.0) + R_pow(ZNKlambdaonefourth, 2.0);
    } else {
      // These are not defined because the c_j's and alpha_j's were not computed
      ZBlambda = R_NaReal;
      ZNKlambdaonefourth = R_NaReal;
      XlambdaAPD = R_NaReal;
    }


    double ZlambdaEPDstar, mutilde, stat = 0.0;
    if (islambdaone) { // lambda = 1
      // Equ. (6.3)
      ZlambdaEPDstar = sqrtn * (Klambda - onemgamma) / sqrt(ctepi);
    } else if (islambdatwo) { // lambda = 2
      // Equ. (6.4)
      ZlambdaEPDstar = sqrtn * (Klambda - ctegamma) / sqrt(ctepibisbis);
    } else {
      // Equ. (6.2)
      ZlambdaEPDstar = sqrtn * (Klambda - ctelambda) / sqrt(ctebeta);
    }
    // Equ. (6.19)
    mutilde = minx;
    double sigmatilde = 0.0;
    for (i = 0; i < n; i++) {
      tmp = R_pow(x[i] - mutilde, lambda);
      sigmatilde = sigmatilde + tmp;
    }
    sigmatilde = R_pow(sigmatilde / (double)n, invlambda);
    double Klambdatilde = 0.0, ZlambdaHEPDstar;
    for (i = 0; i < n; i++) {
      tmp = fabs((x[i] - mutilde) / sigmatilde);
      if (tmp > 0.0000000000000001) Klambdatilde = Klambdatilde + R_pow(tmp, lambda) * log(tmp);
    }
    Klambdatilde = Klambdatilde / (double)n;
    ZlambdaHEPDstar = sqrtn * (Klambdatilde - ctelambda) / sqrt(ctebeta);


    if (version == 1) {
      stat = XlambdaAPD;
    } else if (version == 2) {
      stat = ZNKlambdaonefourth;
    } else if (version == 3) {
      stat = XlambdaAPDstar;
    } else if (version == 4) {
      stat = ZBlambda;
    } else if (version == 5) {
      stat = ZlambdaEPDstar;
    } else { // if (version == 6) {
      stat = ZlambdaHEPDstar;
    }

	
    statistic[0] = stat; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue106.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
      if (usecrit[0] == 1) { // We use the provided critical values
	if (alter[0] == 0) {
	  if (statistic[0] > critvalR[i] || statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0; // two-sided
	} else if (alter[0] == 1) {
	  if (statistic[0] < critvalL[i]) decision[i] = 1; else decision[i] = 0;  // less
	} else if (alter[0] == 2) {
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0;   // greater
	} else if (alter[0] == 3) {
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
	}
      } else {
	if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
      }
    }
    
// If applicable, we free the unused array of pointers
 delete[] info;
 delete[] Tol;
 delete[] Maxit;

}

// We return
    return;
   
        
  }

  int sgn106(double v) {
    return (v > 0) - (v < 0);
  }

  double myf106(double x, void *info) {
    double lambda = (double)((double*)info)[0];
    int i, n = (int)((double*)info)[1];
    double tmp = 0.0;
    for (i = 0; i < n; i++) {
      tmp = tmp + R_pow(fabs((double)((double*)info)[i + 2] - x), lambda - 1.0) * (double)sgn106((double)((double*)info)[i + 2] - x);
    }
    return(tmp);
  }

  double R_zeroin2(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit)				/* Max # of iterations */
{
    double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
    double tol;
    int maxit;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = * Tol;

    /* First test if we have found a root at an endpoint */
    if(fa == 0.0) {
	*Tol = 0.0;
	*Maxit = 0;
	return a;
    }
    if(fb ==  0.0) {
	*Tol = 0.0;
	*Maxit = 0;
	return b;
    }

    while(maxit--)		/* Main iteration loop	*/
    {
	double prev_step = b-a;		/* Distance from the last but one
					   to the last approximation	*/
	double tol_act;			/* Actual tolerance		*/
	double p;			/* Interpolation step is calcu- */
	double q;			/* lated in the form p/q; divi-
					 * sion operations is delayed
					 * until the last moment	*/
	double new_step;		/* Step at this iteration	*/

	if( fabs(fc) < fabs(fb) )
	{				/* Swap data for b to be the	*/
	    a = b;  b = c;  c = a;	/* best approximation		*/
	    fa=fb;  fb=fc;  fc=fa;
	}
	tol_act = 2*EPSILON*fabs(b) + tol/2;
	new_step = (c-b)/2;

	if( fabs(new_step) <= tol_act || fb == (double)0 )
	{
	    *Maxit -= maxit;
	    *Tol = fabs(c-b);
	    return b;			/* Acceptable approx. is found	*/
	}

	/* Decide if the interpolation can be tried	*/
	if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	    && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
					 * Interpolation may be tried	*/
	    double t1,cb,t2;
	    cb = c-b;
	    if( a==c ) {		/* If we have only two distinct	*/
					/* points linear interpolation	*/
		t1 = fb/fa;		/* can only be applied		*/
		p = cb*t1;
		q = 1.0 - t1;
	    }
	    else {			/* Quadric inverse interpolation*/

		q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
		p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
		q = (q-1.0) * (t1-1.0) * (t2-1.0);
	    }
	    if( p>(double)0 )		/* p was calculated with the */
		q = -q;			/* opposite sign; make p positive */
	    else			/* and assign possible minus to	*/
		p = -p;			/* q				*/

	    if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
		&& p < fabs(prev_step*q/2) )	/* and isn't too large	*/
		new_step = p/q;			/* it is accepted
						 * If p/q is too large then the
						 * bisection procedure can
						 * reduce [b,c] range to more
						 * extent */
	}

	if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
	    if( new_step > (double)0 )	/* than tolerance		*/
		new_step = tol_act;
	    else
		new_step = -tol_act;
	}
	a = b;	fa = fb;			/* Save the previous approx. */
	b += new_step;	fb = (*f)(b, info);	/* Do step to a new approxim. */
	if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
	    /* Adjust c for it to have a sign opposite to that of b */
	    c = a;  fc = fa;
	}

    }
    /* failed! */
    *Tol = fabs(c-b);
    *Maxit = -1;
    return b;
}
}
