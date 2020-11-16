#include "aggregator.h"

void Aggregator::add_rows(BendersCut &cut)
{

#pragma omp parallel for
  for (size_t s = 0; s < d_cgmips.size(); ++s)
    d_cgmips[s].add_row(cut);

#pragma omp parallel for
  for (size_t s = 0; s < d_trees.size(); ++s)
    d_trees[s].add_benders_cut(cut);

#pragma omp parallel for
  for (size_t s = 0; s < d_zk.size(); ++s)
    d_zk[s].add_cglp_row(cut);

#pragma omp parallel for
  for (size_t s = 0; s < d_lr.size(); ++s)
    d_lr[s].add_cut(cut);

}