#include "aggregator.h"

BendersCut Aggregator::bac_cut(Master::Solution sol, Master &mp, double tol, size_t maxRounds, double rho_tol)
{
  double rho = sol.thetaVal;
  double *x = sol.xVals.data();
  double cRho = rho_tol + 1;
  BendersCut cut;

  while (cRho > rho_tol)
  {
    cRho = -rho;
    cut = BendersCut{ 0, vector<double>(d_n1), 0};

#pragma omp parallel for reduction(sum : cut) reduction(+:cRho)
    for (size_t s = 0; s < d_trees.size(); ++s)
    {
      double prob = d_probs[s];
      cut += d_trees[s].generate_cut(x, rho, mp, maxRounds, true) * prob;
      cRho -= prob * d_trees[s].cglp_val();
    }


    rho += cRho / (1 + cut.d_tau);
  }

  return cut;
}