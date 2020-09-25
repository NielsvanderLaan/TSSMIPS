#include "pslp.h"

BendersCut Pslp::best_zk_cut(Master::Solution sol, Master &master, size_t maxRounds, bool lap_cuts, double tol)
{
  double *x = sol.xVals.data();
  double rho = sol.thetaVal;
  double cRho = 1;
  
  for (size_t s = 0; s != d_S; ++s)        // generate cutting planes
  {
    if (lap_cuts)
    {
      d_zk[s].update(x, GRB_INFINITY);
      d_zk[s].solve(x, GRB_INFINITY, master, maxRounds, false, false);
    } else
    {
      d_zk[s].update(x, rho);
      d_zk[s].solve(x, rho, master, maxRounds, true, true);
    }
  }

  BendersCut cut;
  
  while (cRho > tol)
  {
    cRho = -rho;
    cut = BendersCut{ 0, vector<double>(d_n1), 0 };

    for (size_t s = 0; s != d_S; ++s)
    {
      d_zk[s].update(x, rho);
      d_zk[s].optimize();
      cut += d_zk[s].subgradient() * d_probs[s];

      cRho += d_probs[s] * d_zk[s].d_objVal;
    }
    
    rho += cRho / (1 + cut.d_tau);
  }
  return cut;
}



















