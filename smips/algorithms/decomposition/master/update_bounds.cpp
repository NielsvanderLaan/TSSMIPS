#include "master.h"

void Master::update_bounds(size_t var, double val, bool lower)
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
}