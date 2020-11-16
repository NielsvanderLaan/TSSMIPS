#include "benders.h"

Benders::Benders(GRBEnv &env, GRBenv *c_env, Problem &problem, vector<Type> &types, bool zk_safe)
:
d_problem(problem),
d_env(env),
d_n1(problem.d_n1),
d_p1(problem.d_p1),
d_n2(problem.d_n2),
d_m2(problem.d_m2),
d_S(problem.d_S),
d_master(env, c_env, problem, zk_safe),
d_gomory(env, problem),
d_agg(env, c_env, problem, types),
d_lb(problem.d_l1),
d_ub(problem.d_u1),
d_visited(problem.d_S),
d_objectives(problem.d_S),
d_UB(GRB_INFINITY)
{
  d_xvals = new double[d_n1];
  d_incumbent = new double[d_n1];
}