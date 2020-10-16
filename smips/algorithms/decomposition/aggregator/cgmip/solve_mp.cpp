#include "cgmip.h"

bool CGMip::solve_mp(bool focus, bool affine, double M)
{
  cout << "solve_mp() called" << endl;
  if (focus)
  {
    d_mp.reset();
    d_mp.set(GRB_IntParam_ScaleFlag, 0);
    d_mp.set(GRB_IntParam_NumericFocus, 3);
    d_mp.set(GRB_IntParam_Method, 1);
    d_mp.set(GRB_IntParam_Presolve, 0);
  }
  cout << "call to optimize()" << endl;
  d_mp.optimize();
  cout << "done" << endl;

  if (mp_optimal())
  {
    if (mp_max_coeff() > M)
    {
      set_mp_bounds(M, affine);
      cout << "call to optimize()" << endl;
      d_mp.optimize();
      cout << "done" << endl;
      set_mp_bounds(GRB_INFINITY, affine);
    }
  }

  if (focus)    // reset params to default
  {
    d_mp.set(GRB_IntParam_ScaleFlag, -1);
    d_mp.set(GRB_IntParam_NumericFocus, 0);
    d_mp.set(GRB_IntParam_Method, -1);
    d_mp.set(GRB_IntParam_Presolve, -1);
  }

  cout << "solve_mp() finished" << endl;
  return mp_optimal();
}