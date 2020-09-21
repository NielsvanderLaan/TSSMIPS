#include "aggregator.h"

BendersCut Aggregator::strong_cut(Master::Solution sol, vector<double> &vx, bool affine, double tol, bool int_feas, double rho_tol)
{
  double rho = sol.thetaVal;
  double *x = sol.xVals.data();
  double cRho = rho_tol + 1;
  BendersCut cut;

  bool first_time = true;
  size_t iter = 0;
  while (cRho > rho_tol)
  {
    ++iter;
    //cout << "iteration: " << iter << '\n';
    cRho = -rho;
    cut = BendersCut{ 0, vector<double>(d_n1), 0};
    double gap = 0;
    for (size_t s = 0; s != d_cgmips.size(); ++s)
    {
      //cout << "s = " << s << '\n';
      double prob = d_probs[s];
      cut += d_cgmips[s].generate_cut(x, rho, first_time, vx[s], affine, tol, int_feas, gap) * prob;
      cRho -= prob * d_cgmips[s].mp_val();
    }
    gap /= d_cgmips.size();
    if (gap > tol)
      cout << "strong_cut(), gap: " << gap << '\n';
    if (affine)
      break;

    rho += cRho / (1 + cut.d_tau);

    first_time = false;
  }
  //cout << "tau: " << cut.d_tau << '\n';
  return cut;
}


