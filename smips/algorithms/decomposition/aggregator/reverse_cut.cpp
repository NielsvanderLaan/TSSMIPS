#include "aggregator.h"

void Aggregator::reverse_cut(double UB)
{

#pragma omp parallel for
  for (size_t s = 0; s < d_cgmips.size(); ++s)
    d_cgmips[s].reverse_cut(UB);

#pragma omp parallel for
  for (size_t s = 0; s < d_trees.size(); ++s)
    d_trees[s].reverse_cut(UB);

#pragma omp parallel for
  for (size_t s = 0; s < d_zk.size(); ++s)
    d_zk[s].reverse_cut(UB);

#pragma omp parallel for
  for (size_t s = 0; s < d_lr.size(); ++s)
    d_lr[s].reverse_cut(UB);

}