#include "benders.h"

Benders::Benders(GRBEnv &env, GRBenv *c_env, Problem &problem, bool zk_safe)
:
d_problem(problem),
d_env(env),
d_master(env, c_env, problem, zk_safe),
d_sub(env, problem, -1),
d_lr(env, problem),
d_gomory(env, problem),
d_ald(c_env, problem),
d_pslp(env, c_env, problem),
d_agg(env, c_env, problem),
d_lb(problem.d_l1),
d_ub(problem.d_u1),
d_S(problem.d_S),
d_n1(problem.d_n1),
d_p1(problem.d_p1),
d_n2(problem.d_n2),
d_m2(problem.d_m2),
d_visited(problem.d_S),
d_objectives(problem.d_S)
{
  d_xvals = new double[d_n1];
  d_incumbent = new double[d_n1];
  if (not problem.d_fix_rec)
  {
    for (size_t s = 0; s != d_S; ++s)
      d_sub_omega.push_back(Sub {env, problem, s});
  }
}