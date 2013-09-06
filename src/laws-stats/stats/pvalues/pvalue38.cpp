if (alter[0] == 0) pvalue[0] = 2.0*pnorm5(fabs(statvalue),0.0,1.0,0,0);
if (alter[0] == 1) pvalue[0] = pnorm5(statvalue,0.0,1.0,1,0);
if (alter[0] == 2) pvalue[0] = pnorm5(statvalue,0.0,1.0,0,0); 
