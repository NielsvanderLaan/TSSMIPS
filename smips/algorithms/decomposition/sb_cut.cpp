#include "benders.h" 

#include <omp.h>

BendersCut Benders::sb_cut(double *x)
{
  BendersCut cut {0.0, vector<double> (d_n1, 0.0), 0.0};

  vector<Sub> subs(d_S, d_sub);
  //vector<Lagrangian> lrs(d_S, d_lr);

#pragma omp parallel for reduction(sum : cut) firstprivate(d_lr) num_threads(4)
  for (size_t s = 0; s < d_S; ++s)
  {
    //vector<double> pi = d_sub.compute_slope(s, x);
    vector<double> pi = subs[s].compute_slope(s, x);

    cut += d_lr.strong_cut(s, pi) * d_problem.d_probs[s];
    //cut += lrs[s].strong_cut(s, pi) * d_problem.d_probs[s];

  }

  return cut;
}