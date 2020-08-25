#include "tree.h"

bool Tree::add_branch(size_t node_idx, Split split)
{

  Benders *child = new Benders(*d_nodes[node_idx]);      // left child
  child->update_bounds(split.var, split.left, false);
  d_nodes.push_back(child);
  d_LB_nodes.push_back(d_LB_nodes[node_idx]);   // child inherits LB from parent
  d_nodes[node_idx]->update_bounds(split.var, split.right, true);       // right child

  return false;
}