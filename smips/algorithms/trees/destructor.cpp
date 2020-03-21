#include "tree.h"

Tree::~Tree()
{
  for (Benders *node_ptr : d_nodes)
    delete node_ptr;
}