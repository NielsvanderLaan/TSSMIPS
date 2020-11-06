#include "aggregator.h"

void Aggregator::reverse_cut(double UB)
{
  for (size_t s = 0; s != d_cgmips.size(); ++s)
  {
    d_cgmips[s].reverse_cut(UB);
    d_trees[s].reverse_cut(UB);
    d_zk[s].reverse_cut(UB);
    d_lr[s].reverse_cut(UB);
  }
}