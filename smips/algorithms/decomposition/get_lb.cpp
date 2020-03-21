#include "benders.h"

double Benders::get_lb()
{
  double LB;
  GRBgetdblattr(d_master.d_cmodel, "ObjVal", &LB);  
  return LB;
}