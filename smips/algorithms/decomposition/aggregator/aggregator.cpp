#include "aggregator.h"

Aggregator::Aggregator(GRBEnv &env, GRBenv *c_env, Problem &problem)
:  
  d_n1(problem.d_n1),
  d_probs(problem.d_probs),
  d_vw(env),
  d_problem(problem),
  d_fix_rec(problem.d_fix_rec)
{
  init_vw(problem);

  size_t S = problem.d_S;
  d_cgmips.reserve(S);
  d_trees.reserve(S);
  d_zk.reserve(S);
  d_sub.reserve(S);
  d_lr.reserve(S);

  for (size_t s = 0; s != problem.d_S; ++s)
  {
    /*
    GRBEnv subenv;
    subenv.set(GRB_IntParam_OutputFlag, 0);
    subenv.set(GRB_IntParam_Threads, 1);
    */

    CGMip cgmip{env, problem, s };
    d_cgmips.push_back(cgmip);

    ZK zk{c_env, env, problem, s};
    d_zk.push_back(zk);

    ZkTree tree{ c_env, env, problem, s };
    d_trees.push_back(tree);

    Sub sub {env, problem, s};
    d_sub.push_back(sub);

    Lagrangian lr(env, problem, s);
    d_lr.push_back(lr);
  }
}