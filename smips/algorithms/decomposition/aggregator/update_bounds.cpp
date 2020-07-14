#include "aggregator.h"

void Aggregator::update_bounds(size_t var, double val, bool lower)
{
  for (size_t s = 0; s != d_cgmips.size(); ++s)
  {
    d_cgmips[s].update_bound(var, val, lower);
    d_trees[s].update_fs_bounds(var, val, lower);
  }


}