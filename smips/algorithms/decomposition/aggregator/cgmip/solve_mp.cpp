#include "cgmip.h"

bool CGMip::solve_mp(bool focus)
{
  if (focus)
  {
    d_mp.reset();
    d_mp.set(GRB_IntParam_ScaleFlag, 0);
    d_mp.set(GRB_IntParam_NumericFocus, 3);
  }
  d_mp.optimize();
  if (focus)    // reset params to default
  {
    d_mp.set(GRB_IntParam_ScaleFlag, -1);
    d_mp.set(GRB_IntParam_NumericFocus, 0);
  }
  return mp_optimal();
}