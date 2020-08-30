#include "benders.h" 

BendersCut Benders::sb_cut(double *x)
{
  BendersCut cut {0.0, vector<double> (d_n1, 0.0), 0.0};

  for (size_t s = 0; s != d_S; ++s)
  {
    vector<double> pi = d_sub.compute_slope(s, x);

    cut += d_lr.strong_cut(s, pi) * d_problem.d_probs[s];
  }

  return cut;
}