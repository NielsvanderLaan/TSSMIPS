#include "zktree.h"

ZkTree::~ZkTree()
{
  for (ZK *node : d_nodes)
    delete node;
}