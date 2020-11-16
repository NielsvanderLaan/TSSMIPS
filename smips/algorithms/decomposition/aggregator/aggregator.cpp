#include "aggregator.h"

Aggregator::Aggregator(GRBEnv &env, GRBenv *c_env, Problem &problem, vector<Type> &types)
:  
  d_n1(problem.d_n1),
  d_probs(problem.d_probs),
  d_vw(env),
  d_problem(problem),
  d_fix_rec(problem.d_fix_rec)
{
  init_vw(problem);

  size_t S = problem.d_S;
  if (find(types.begin(), types.end(), SC_RG) != types.end() || find(types.begin(), types.end(), LR) != types.end())
  {
    d_cgmips.reserve(S);
    for (size_t s = 0; s != problem.d_S; ++s)
    {
      CGMip cgmip{env, problem, s };
      d_cgmips.push_back(cgmip);
    }
  }

  if (find(types.begin(), types.end(), SC_BAB) != types.end() || find(types.begin(), types.end(), SC_BAC) != types.end())
  {
    d_trees.reserve(S);
    for (size_t s = 0; s != problem.d_S; ++s)
    {
      ZkTree tree{c_env, env, problem, s};
      d_trees.push_back(tree);
    }
  }

  if (find(types.begin(), types.end(), LR_LAP) != types.end() ||
      find(types.begin(), types.end(), SC_ZK) != types.end()  ||
      find(types.begin(), types.end(), SC_LAP) != types.end() )
  {
    d_zk.reserve(S);
    for (size_t s = 0; s != problem.d_S; ++s)
    {
      ZK zk{c_env, env, problem, s};
      d_zk.push_back(zk);
    }
  }

  if (find(types.begin(), types.end(), SB) != types.end())
  {
    d_lr.reserve(S);
    for (size_t s = 0; s != problem.d_S; ++s)
    {
      Lagrangian lr(env, problem, s);
      d_lr.push_back(lr);
    }
  }

  d_sub.reserve(S);

  for (size_t s = 0; s != problem.d_S; ++s)
  {
    Sub sub {env, problem, s};
    d_sub.push_back(sub);
  }
}