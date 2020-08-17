#include "zktree.h"

bool ZkTree::is_feasible(ZK *node)
{
  return all_of(node->d_yvals.begin(), node->d_yvals.begin() + node->d_p2, [](double val){ return is_integer(val); });
}