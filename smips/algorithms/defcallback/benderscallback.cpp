#include "benderscallback.h"

BendersCallback::BendersCallback(Problem &problem, GRBEnv &env, GRBenv *c_env, DEF &def)
:
  d_xvars(def.d_xvars),
  d_theta(def.d_theta),
  d_ben(env, c_env, problem),
  d_problem(problem),
  d_ncuts(0)
{
  d_ben.lpSolve();
};