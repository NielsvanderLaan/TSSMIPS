#include "aggregator.h"

BendersCut Aggregator::sb_cut(vector<double> &x)
{
  BendersCut cut {0.0, vector<double> (d_n1, 0.0), 0.0};

  BendersCut cuts[d_cgmips.size()];

#pragma omp parallel for reduction(sum : cut)
  for (size_t s = 0; s < d_sub.size(); ++s)
  {
    vector<double> pi = d_sub[s].compute_slope(x);
    cuts[s] = d_lr[s].strong_cut(pi) * d_problem.d_probs[s];
    //cut += d_lr[s].strong_cut(pi) * d_problem.d_probs[s];
  }

  for (size_t s = 0; s < d_cgmips.size(); ++s)
    cut += cuts[s];


  return cut;
}
