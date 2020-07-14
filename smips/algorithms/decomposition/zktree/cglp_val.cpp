#include "ZkTree.h"

double ZkTree::cglp_val()
{
  if (d_cglp.get(GRB_IntAttr_Status) == 2)
    return d_cglp.get(GRB_DoubleAttr_ObjVal);

  return GRB_INFINITY;
}