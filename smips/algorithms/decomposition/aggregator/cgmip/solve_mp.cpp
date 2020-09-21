#include "cgmip.h"

bool CGMip::solve_mp(bool focus, double M)
{
  if (focus)
  {
    d_mp.reset();
    d_mp.set(GRB_IntParam_ScaleFlag, 0);
    d_mp.set(GRB_IntParam_NumericFocus, 3);
    d_mp.set(GRB_IntParam_Method, 0);
  }

  d_mp.optimize();

  if (mp_optimal())
  {
    if (mp_max_coeff() > M)
    {
      set_mp_bounds(M);
      d_mp.optimize();
      set_mp_bounds(GRB_INFINITY);
    }
  }

  if (focus)    // reset params to default
  {
    d_mp.set(GRB_IntParam_ScaleFlag, -1);
    d_mp.set(GRB_IntParam_NumericFocus, 0);
    d_mp.set(GRB_IntParam_Method, -1);
  }

  return mp_optimal();
}