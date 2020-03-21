#include "pslp.h"

Pslp::Pslp(GRBEnv &env, GRBenv *c_env, Problem &problem)
:
d_n1(problem.d_n1),
d_S(problem.d_S),
d_probs(problem.d_probs)
{
  d_zk.reserve(d_S);
  
  for (size_t scenario = 0; scenario != d_S; ++scenario)
    d_zk.emplace_back(ZK{c_env, env, problem, scenario});
}