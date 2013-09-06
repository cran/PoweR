pval = exp(-7.01256 * Kd*Kd * ((double)nd + 2.78019) + 2.99587 * Kd * sqrt((double)nd + 2.78019) - 0.122119 + 0.974598/sqrt((double)nd) + 1.67997/(double)nd);

if (pval > 0.1) {
  KK = (sqrt((double)n) - 0.01 + 0.85/sqrt((double)n)) * statKS;
  if (KK <= 0.302) {
    pval = 1.0;
  } else if (KK <= 0.5) {
    pval = 2.76773 - 19.828315 * KK + 80.709644 * KK*KK - 138.55152 * KK*KK*KK + 81.218052 * KK*KK*KK*KK;
  } else if (KK <= 0.9) {
    pval = -4.901232 + 40.662806 * KK - 97.490286 * KK*KK + 94.029866 * KK*KK*KK - 32.355711 * KK*KK*KK*KK;
  } else if (KK <= 1.31) {
    pval = 6.198765 - 19.558097 * KK + 23.186922 * KK*KK - 12.234627 * KK*KK*KK + 2.423045 * KK*KK*KK*KK;
  } else {
    pval = 0.0;
  }
}

pvalue[0] = pval;
