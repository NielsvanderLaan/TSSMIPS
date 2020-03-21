#include "zktree.h"

ZkTree::ZkTree(GRBenv *env, GRBEnv &cpp_env, Problem &problem, size_t scenario)
{
  ZK *root = new ZK{env, cpp_env, problem, scenario};
  d_nodes.push_back(root);
  
  // TODO
  // initialize CGLP (and d_lambda, d_constrs)
}