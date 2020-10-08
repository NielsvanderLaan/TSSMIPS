#include "master.h"

void Master::chg_mp_tol(bool focus)
{
  if (focus) GRBreset(d_cmodel, 0);

  GRBenv *env = GRBgetenv(d_cmodel);
  GRBsetintparam(env, "NumericFocus", focus ? 3 : 0);
  GRBsetintparam(env, "ScaleFlag", focus ? 0 : -1);
  GRBsetintparam(env, "Method", focus ? 1 : -1);
  GRBsetintparam(env, "Presolve", focus ? 0 : -1);

  GRBupdatemodel(d_cmodel);
}