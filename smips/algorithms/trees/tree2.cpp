#include "tree.h"

Tree::Tree(GRBEnv &env, GRBenv *c_env, Problem &problem, Benders &root)
:
  d_problem(problem),
  d_c_env(c_env),
  d_env(env),
  d_UB_global(1e20)
{
  Benders *root_ptr = new Benders(root);
  d_nodes.push_back(root_ptr);
  double lb = root.get_lb();
  d_LB_nodes.push_back(lb);
  d_LB_global = lb;
}