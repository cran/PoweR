if (alter[0] == 0) {
  statT = 0.0 - fabs(statT);
  pvaltmp = 2.0 * pt(statT, n-1, 1, 0);
 }
 else if (alter[0] == 1) {
   pvaltmp = pt(statT, n-1, 1, 0);
 }
 else if (alter[0] == 2) {
   pvaltmp = pt(statT, n-1, 0, 0);
 }

pvalue[0] = pvaltmp;
