#include "benders.h"

void Benders::update_bounds(size_t var, double val, bool lower)
{
  if (lower)
    d_lb[var] = val;
  else
    d_ub[var] = val;
  
  d_master.update_bounds(var, val, lower);
  d_pslp.update_bounds(var, val, lower);
  d_agg.update_bounds(var, val, lower);
}