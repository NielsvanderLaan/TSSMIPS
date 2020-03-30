#include "benders.h"

double Benders::get_lb()
{
  double LB;
  GRBgetdblattr(d_master.d_cmodel, "ObjBound", &LB);
  return LB;
}