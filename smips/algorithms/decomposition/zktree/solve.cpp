#include "zktree.h"

bool ZkTree::solve(size_t node_idx, double *x, double theta, Master &master, size_t maxRounds, bool gomory)
{

  //bool feasible = d_nodes[node_idx]->solve(x, theta, master, maxRounds, gomory);    // BAC

  bool feasible = d_nodes[node_idx]->optimize(); // pure BAB

  // TODO
  // update cglp (include any new cuts)
  
  return feasible;
} 