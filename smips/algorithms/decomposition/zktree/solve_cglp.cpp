#include "zktree.h"

void ZkTree::solve_cglp(double M)
{
  optimize();

  if (max_coeff(d_candidate) > 1e8)
  {
    cglp_bounds(M, true);
    optimize();
    cglp_bounds(GRB_INFINITY, false);
  }

}