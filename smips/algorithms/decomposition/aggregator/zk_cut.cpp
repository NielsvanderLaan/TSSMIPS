#include "aggregator.h"

BendersCut Aggregator::zk_cut(Master::Solution sol, Master &master, bool lap_cuts, bool affine, size_t maxRounds, double tol, double rho_tol)
{
  double *x = sol.xVals.data();
  double theta = sol.thetaVal;
  double rho = theta;
  double cRho = tol + 1;    // this choice ensures that the while loop is always entered

  BendersCut cut;

  while (cRho > rho_tol)
  {
    cRho = -rho;
    cut = BendersCut{ 0, vector<double>(d_n1, 0.0), 0 };

#pragma omp parallel for reduction(sum : cut) reduction(+:cRho)
    for (size_t s = 0; s < d_zk.size(); ++s)
    {
      d_zk[s].update(x,affine ? GRB_INFINITY : rho);
      d_zk[s].solve(x, theta, affine ? GRB_INFINITY : rho, master, maxRounds, not lap_cuts);

      cut += d_zk[s].subgradient() * d_probs[s];
      cRho += d_probs[s] * d_zk[s].d_objVal;
    }

    if (affine)
      break;

    rho += cRho / (1 + cut.d_tau);
  }

#pragma omp parallel for
  for (size_t s = 0; s < d_zk.size(); ++s)
    d_zk[s].clear();      // reduces memory usage


  return cut;
}

