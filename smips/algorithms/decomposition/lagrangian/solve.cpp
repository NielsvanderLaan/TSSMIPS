#include "lagrangian.h"

double Lagrangian::solve()
{
  d_model.optimize();
  if (d_model.get(GRB_IntAttr_Status) != 2)
  {
    d_model.write("lr.lp");
    exit(159);
  }
  return d_model.get(GRB_DoubleAttr_ObjVal);
}