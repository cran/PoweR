double *xoptim, *loweroptim, *upperoptim, *Fminoptim, *exoptim;
int *nbdoptim, *failoptim, *fncountoptim, *grcountoptim;
char *msgoptim;

xoptim = new double[1];
loweroptim  = new double[1];
upperoptim = new double[1];
Fminoptim = new double[1];
exoptim = new double[2];

nbdoptim = new int[1];
failoptim = new int[1];
fncountoptim = new int[1];
grcountoptim = new int[1];

msgoptim = new char[60]; // Always 60.

xoptim[0] = 0.3819660112501051; // initial value (a + (1-phi)(b-a) where ‘(a,b) = (lower, upper)’ and phi = (sqrt(5) - 1)/2 = 0.61803..  is the golden section ratio.), and solution found on exit
loweroptim[0] = 0.00000000000000000001; // lower bound
upperoptim[0] = 1.0; // upper bound
nbdoptim[0] = 2; // 0 if x(i) is unbounded, 1 if x(i) has only a lower bound, 2 if x(i) has both lower and upper bounds, 3 if x(i) has only an upper bound.
Fminoptim[0] = 0.0; // To display some information during the calculations
failoptim[0] = 0; // Contains some value indicating what kind of error occured
exoptim[0] = (double)n; // Value to be passed to the function quantileestimC
exoptim[1] = statvalue; // Value to be passed to the function quantileestimC
double factroptim = 100000000.0; // decrease for higher accuracy. The iteration will stop when (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch where epsmch is the machine precision. Typical values for factr: 1.d+12 for low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely high accuracy.
double pgtoloptim = 0.0; // The iteration will stop when max{|proj g_i | i = 1, ..., n} <= pgtol where pg_i is the ith component of the projected gradient.
fncountoptim[0] = 0; 
grcountoptim[0] = 0;
int maxitoptim = 100;
msgoptim[0] = '\0'; // Used to display some messages durint the calculations.
int traceoptim = 0; // To display some information during the calculations
/*
If 2: iprint = 0
If 3: iprint = nREPORT (voir plus bas)
If 4: iprint = 99
If 5: iprint = 100
If 6: iprint = 101
Sinon: iprint = -1

       iprint is an integer variable that must be set by the user.
	 It controls the frequency and type of output generated:
	  iprint<0    no output is generated;
	  iprint=0    print only one line at the last iteration;
	  0<iprint<99 print also f and |proj g| every iprint iterations;
	  iprint=99   print details of every iteration except n-vectors;
	  iprint=100  print also the changes of active set and final x;
	  iprint>100  print details of every iteration including x and g;
	 When iprint > 0, the file iterate.dat will be created to
			  summarize the iteration.
 */
int nREPORToptim = 10;
/*
 */

//        This subroutine solves bound constrained optimization problems by
//	 using the compact formula of the limited memory BFGS updates.
// First 1 is for number of variables
// Second 5 is for the maximum number of variable metric corrections allowed in the limited memory matrix.
lbfgsb(1, 5, xoptim, loweroptim,
	    upperoptim, nbdoptim, Fminoptim, quantileestimC,
	    myoptimgrC, failoptim, exoptim, factroptim,
	    pgtoloptim, fncountoptim, grcountoptim,
	    maxitoptim, msgoptim, traceoptim, nREPORToptim);

pvalue[0] = xoptim[0];

double M;
double *xtmp;
xtmp = new double[n];
// calculate sample median
for (i=0;i<=(n-1);i++) {
  xtmp[i] = x[i];
 }	  
R_rsort(xtmp,n); // We sort the data
if ((n%2) == 0) {
  M = (xtmp[n/2]+xtmp[n/2-1])/2.0; 
 } else {
  M = xtmp[n/2]; // sample median
 }
delete[] xtmp;


if (alter[0] == 0) {if (statvalue < M) pvalue[0] = 2.0 * pvalue[0]; else pvalue[0] = 2.0 * (1.0 - pvalue[0]);}
if (alter[0] == 2) pvalue[0] = 1.0 - pvalue[0]; 
/*
// Avant, j'avais mis cela:
// if (alter[0] == 0) pvalue[0] = 2.0 * pnorm5(fabs(statvalue), 0.0, 1.0, 0, 0);
// if (alter[0] == 1) pvalue[0] = pnorm5(statvalue, 0.0, 1.0, 1, 0);
// if (alter[0] == 2) pvalue[0] = pnorm5(statvalue, 0.0, 1.0, 0, 0); 
*/

delete[] xoptim;
delete[] loweroptim;
delete[] upperoptim;
delete[] Fminoptim;
delete[] exoptim;

delete[] nbdoptim;
delete[] failoptim;
delete[] fncountoptim;
delete[] grcountoptim;

delete[] msgoptim;





