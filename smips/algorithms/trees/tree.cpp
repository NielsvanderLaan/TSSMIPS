#include "tree.h"

Tree::Tree(GRBEnv &env, GRBenv *c_env, Problem &problem)
:
  d_problem(problem),
  d_c_env(c_env),
  d_env(env),
  d_UB_global(GRB_INFINITY)
{
  Benders *root = new Benders(env, c_env, problem, false);
  double lb = root->lpSolve();
  d_nodes.push_back(root);
  d_LB_nodes.push_back(lb);
  d_LB_global = lb;
}