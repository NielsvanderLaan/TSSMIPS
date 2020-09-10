#include "aggregator.h"

BendersCut Aggregator::strong_cut(Master::Solution sol, vector<double> &vx, bool affine, double tol, double rho_tol)
{
  double rho = sol.thetaVal;
  double *x = sol.xVals.data();
  double cRho = rho_tol + 1;
  BendersCut cut;

  bool first_time = true;
  while (cRho > rho_tol)
  {
    cRho = -rho;
    cut = BendersCut{ 0, vector<double>(d_n1), 0};
    double gap = 0;
    for (size_t s = 0; s != d_cgmips.size(); ++s)
    {
      double prob = d_probs[s];
      cut += d_cgmips[s].generate_cut(x, rho, first_time, vx[s], affine, tol, gap) * prob;
      cRho -= prob * d_cgmips[s].mp_val();
    }
    cout << "strong_cut(), gap: " << gap / d_cgmips.size() << '\n';
    if (affine)
      break;

    rho += cRho / (1 + cut.d_tau);

    first_time = false;
  }
  //cout << "tau: " << cut.d_tau << '\n';
  return cut;
}


