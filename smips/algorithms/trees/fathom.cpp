#include "tree.h"

void Tree::fathom()
{
  for (size_t node_idx = d_nodes.size() - 1; node_idx != -1; --node_idx)
  {
    if (d_LB_nodes[node_idx] > d_UB_global + 1e-8)    // fathom node (safe side here)
    {
      delete d_nodes[node_idx];
      d_nodes.erase(d_nodes.begin() + node_idx);
      d_LB_nodes.erase(d_LB_nodes.begin() + node_idx);
    }
  }
}