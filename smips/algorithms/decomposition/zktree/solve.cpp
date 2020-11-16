#include "zktree.h"

bool ZkTree::solve(size_t node_idx, double *x, double theta, double rho, Master &master, bool cuts, size_t maxRounds, double tol)
{
  if (cuts)
    return d_nodes[node_idx]->solve(x, theta, rho, master, maxRounds, true, tol);    // BAC
  else
    return d_nodes[node_idx]->optimize();                                                 // pure BAB
} 