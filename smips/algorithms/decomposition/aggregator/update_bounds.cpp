#include "aggregator.h"

void Aggregator::update_bounds(size_t var, double val, bool lower)
{
#pragma omp parallel for
  for (size_t s = 0; s < d_cgmips.size(); ++s)
    d_cgmips[s].update_bound(var, val, lower);

#pragma omp parallel for
  for (size_t s = 0; s < d_trees.size(); ++s)
    d_trees[s].update_fs_bounds(var, val, lower);

#pragma omp parallel for
  for (size_t s = 0; s < d_zk.size(); ++s)
    d_zk[s].update_bound(var, val, lower, true);

#pragma omp parallel for
  for (size_t s = 0; s < d_lr.size(); ++s)
    d_lr[s].update_bound(var, val, lower);

}