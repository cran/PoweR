if (alter[0] == 0) pvalue[0] = 2.0 * Rf_pnorm5(fabs(stat84),0.0,1.0,0,0);
if (alter[0] == 1) pvalue[0] = Rf_pnorm5(stat84,0.0,1.0,1,0);
if (alter[0] == 2) pvalue[0] = Rf_pnorm5(stat84,0.0,1.0,0,0);

