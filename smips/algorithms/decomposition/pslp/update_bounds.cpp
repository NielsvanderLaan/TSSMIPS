#include "pslp.h"

void Pslp::update_bounds(size_t var, double val, bool lower)
{
  for (size_t s = 0; s != d_S; ++s)
    d_zk[s].update_bound(var, val, lower, true);
}