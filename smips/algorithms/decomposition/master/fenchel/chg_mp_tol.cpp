#include "fenchel.h"

void Fenchel::chg_mp_tol(bool focus)
{
  if (focus)  d_mp.reset();

  d_mp.set(GRB_IntParam_ScaleFlag, focus ? 0 : -1);
  d_mp.set(GRB_IntParam_NumericFocus, focus ? 3 : 0);
  d_mp.set(GRB_IntParam_Method, focus ? 1 : -1);
  d_mp.set(GRB_IntParam_Presolve, focus ? 0 : -1);
}