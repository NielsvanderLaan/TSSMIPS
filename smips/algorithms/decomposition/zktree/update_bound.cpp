#include "zktree.h"

void ZkTree::update_bound(size_t node_idx, size_t var_idx, bool lower)
{
  ZK *node = d_nodes[node_idx];
  vector<double> &yvals = d_nodes[node_idx]->d_yvals;
  double val = lower ? ceil(yvals[var_idx]) : floor(yvals[var_idx]);
  node->update_bound(var_idx, val, lower, false);
  
  // TODO
  // update cglp
}