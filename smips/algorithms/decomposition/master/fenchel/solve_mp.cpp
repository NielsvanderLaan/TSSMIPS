#include "fenchel.h"

bool Fenchel::solve_mp()
{
  d_mp.optimize();

  return d_mp.get(GRB_IntAttr_Status) == 2;
}