#include "aggregator.h"

BendersCut Aggregator::strong_cut(Master::Solution sol, vector<double> &vx, bool affine, double tol, double rho_tol)
{
  cout << "strong_cut()\n";
  double rho = sol.thetaVal;
  double *x = sol.xVals.data();
  double cRho = rho_tol + 1;
  BendersCut cut;

  bool first_time = true;
  while (cRho > rho_tol)
  {
    cout << "rho = " << rho << '\n';
    cRho = -rho;
    cut = BendersCut{ 0, vector<double>(d_n1), 0};

    for (size_t s = 0; s != d_cgmips.size(); ++s)
    {
      cout << "s = " << s << '\n';
      double prob = d_probs[s];
      cut += d_cgmips[s].generate_cut(x, rho, first_time, vx[s], affine, tol) * prob;
      cRho -= prob * d_cgmips[s].mp_val();
    }

    rho += cRho / (1 + cut.d_tau);

    first_time = false;
  }

  //cout << "tau: " << cut.d_tau << '\n';
  return cut;
}


