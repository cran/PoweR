pvalue[0] = 0;
if (version == 1) {
  pvalue[0] = Rf_pchisq(statistic[0], 2.0, 0, 0);
 } else if (version == 2) {
  if (alter[0] == 0) { // 0: two.sided=bilateral
    pvalue[0] = 2.0 * Rf_pnorm5(fabs(statistic[0]), 0.0, 1.0, 0, 0);
  }
  if (alter[0] == 1) { // 1: less=unilateral
    pvalue[0] = Rf_pnorm5(-fabs(statistic[0]), 0.0, 1.0, 1, 0);
  }
  if (alter[0] == 2) { // 2: greater=unilateral
    pvalue[0] = Rf_pnorm5(fabs(statistic[0]), 0.0, 1.0, 0, 0);
  }
 } else if (version == 3) {
  if (alter[0] == 0) { // 0: two.sided=bilateral
    pvalue[0] = 2.0 * Rf_pnorm5(fabs(statistic[0]), 0.0, 1.0, 0, 0);
  }
  if (alter[0] == 1) { // 1: less=unilateral
    pvalue[0] = Rf_pnorm5(-fabs(statistic[0]), 0.0, 1.0, 1, 0);
  }
  if (alter[0] == 2) { // 2: greater=unilateral
    pvalue[0] = Rf_pnorm5(fabs(statistic[0]), 0.0, 1.0, 0, 0);
  }
 } else if (version == 4) {
  pvalue[0] = Rf_pchisq(statistic[0], 2.0, 0, 0);
 } else if (version == 5) {
  if (alter[0] == 0) { // 0: two.sided=bilateral
    pvalue[0] = 2.0 * Rf_pnorm5(fabs(statistic[0]), 0.0, 1.0, 0, 0);
  }
  if (alter[0] == 1) { // 1: less=unilateral
    pvalue[0] = Rf_pnorm5(-fabs(statistic[0]), 0.0, 1.0, 1, 0);
  }
  if (alter[0] == 2) { // 2: greater=unilateral
    pvalue[0] = Rf_pnorm5(fabs(statistic[0]), 0.0, 1.0, 0, 0);
  }
 } else if (version == 6) {
  if (alter[0] == 0) { // 0: two.sided=bilateral
    pvalue[0] = 2.0 * Rf_pnorm5(fabs(statistic[0]), 0.0, 1.0, 0, 0);
  }
  if (alter[0] == 1) { // 1: less=unilateral
    pvalue[0] = Rf_pnorm5(-fabs(statistic[0]), 0.0, 1.0, 1, 0);
  }
  if (alter[0] == 2) { // 2: greater=unilateral
    pvalue[0] = Rf_pnorm5(fabs(statistic[0]), 0.0, 1.0, 0, 0);
  }
 }
