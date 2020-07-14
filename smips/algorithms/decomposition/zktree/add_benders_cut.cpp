#include "zktree.h"

void ZkTree::add_benders_cut(const BendersCut &cut)
{
  for (size_t node_idx = 0; node_idx != d_nodes.size(); ++node_idx)
    add_row_to_cglp(cut.d_beta.data(), 1 + cut.d_tau, 0.0, cut.d_alpha, node_idx);

  d_cglp.update();
}