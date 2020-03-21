#include "zktree.h"

void ZkTree::add_child(size_t node_idx)
{
  d_nodes.emplace_back(new ZK(*d_nodes[node_idx]));
  
  // TODO
  // update cglp
}