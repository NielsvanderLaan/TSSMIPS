#include "master.h"

void Master::update_bounds(size_t var, double val, bool lower)
{
  if (d_zk_safe)
  {
    if (lower)
    {
      GRBsetdblattrelement(d_cmodel, "RHS", var, val);
      d_gamma[var] = val;
    } else
    {
      var += d_n1;
      GRBsetdblattrelement(d_cmodel, "RHS", var, val);
      d_gamma[var] = -val;
    }
  } else  // not zk-safe
  {
    ++var;
    if (lower)
      GRBsetdblattrelement(d_cmodel, "LB", var, val);
    else
      GRBsetdblattrelement(d_cmodel, "UB", var, val);
  }
  GRBupdatemodel(d_cmodel);
}