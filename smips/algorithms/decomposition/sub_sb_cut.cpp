#include "benders.h"

BendersCut Benders::sub_sb_cut(size_t s, double *x, Sub &sub, Lagrangian &lr)
{
  vector<double> pi = sub.compute_slope(s, x);

  return lr.strong_cut(s, pi);
}