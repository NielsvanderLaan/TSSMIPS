#include "zktree.h"

void ZkTree::add_benders_cut(double *coeff_x, double coeff_theta, double rhs)
{
  for (size_t node_idx = 0; node_idx != d_nodes.size(); ++node_idx)
    add_row_to_cglp(coeff_x, coeff_theta, 0.0, rhs, node_idx);
}