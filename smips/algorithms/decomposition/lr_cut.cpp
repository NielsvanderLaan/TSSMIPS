#include "benders.h"

BendersCut Benders::lr_cut(double *x, vector<double> &vx)
{
  BendersCut cut {0.0, vector<double> (d_n1, 0.0), 0.0};

  for (size_t s = 0; s != d_S; ++s)
  {
    vector<double> pi = d_sub.compute_slope(s, x);
    pi = vector<double> (pi.size(), 0.0);
    cut += d_lr.lr_cut(s, x, vx[s], pi) * d_problem.d_probs[s];
  }

  return cut;
}