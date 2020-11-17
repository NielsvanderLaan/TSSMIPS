#include "zktree.h"

void ZkTree::clear()
{
  for (ZK *node : d_nodes)
    node->clear();
}