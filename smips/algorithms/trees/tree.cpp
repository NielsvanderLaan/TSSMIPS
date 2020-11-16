#include "tree.h"

Tree::Tree(GRBEnv &env, GRBenv *c_env, Problem &problem, vector<Type> &types)
:
  d_problem(problem),
  d_c_env(c_env),
  d_env(env),
  d_UB_global(GRB_INFINITY)
{
  Benders *root = new Benders(env, c_env, problem, types, true);

  double lb = root->lpSolve();
  //double lb = root->strong_benders();
  d_nodes.push_back(root);
  d_LB_nodes.push_back(lb);
  d_LB_global = lb;
}