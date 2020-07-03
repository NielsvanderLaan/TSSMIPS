#include "zktree.h"

ZkTree::ZkTree(GRBenv *env, GRBEnv &cpp_env, Problem &problem, size_t scenario)
:
  d_cglp(cpp_env),
  d_L(problem.d_L),
  d_beta(problem.d_n1),
  d_constrs(1),
  d_lambda(1),
  d_lb_mult_inds{vector<int>(problem.d_n1, -1)},
  d_ub_mult_inds{vector<int>(problem.d_n1, -1)}
{

  ZK *root = new ZK{env, cpp_env, problem, scenario};
  d_nodes.push_back(root);

  // initializing cglp
  d_alpha = d_cglp.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

  vector<double> beta_lb(problem.d_n1, -GRB_INFINITY);
  GRBVar *beta_ptr = d_cglp.addVars(beta_lb.data(), NULL, NULL, NULL, NULL, problem.d_n1);
  copy_n(beta_ptr, problem.d_n1, d_beta.begin());
  delete[] beta_ptr;

  d_tau = d_cglp.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
  d_kappa = d_cglp.addVar(1.0, 1.0, 0.0, GRB_CONTINUOUS);

  size_t m1 = problem.d_m1;
  vector<double> lambda_lb(m1, -GRB_INFINITY);
  vector<double> lambda_ub(m1, GRB_INFINITY);
  fill_n(lambda_lb.begin() + problem.d_fs_leq, problem.d_fs_geq, 0.0);
  fill_n(lambda_ub.begin(), problem.d_fs_leq, 0.0);

  GRBVar *lambda = d_cglp.addVars(lambda_lb.data(), lambda_ub.data(), NULL, NULL, NULL, m1);
  for (size_t mult = 0; mult != m1; ++mult)
    d_lambda[0].push_back(lambda[mult]);

  GRBLinExpr alpha_lhs;
  alpha_lhs.addTerms(problem.d_b.data(), lambda, m1);

  for (size_t var = 0; var != d_beta.size(); ++var)
  {
    double col[m1];
    for (size_t row = 0; row != m1; ++row)
      col[row] = problem.d_Amat[row][var];
    GRBLinExpr lhs;
    lhs.addTerms(col, lambda, m1);

    if (problem.d_u1[var] < 1e10)
    {
      GRBVar mult = d_cglp.addVar(-GRB_INFINITY, 0.0, 0.0, GRB_CONTINUOUS);
      d_lambda[0].push_back(mult);
      d_ub_mult_inds[0][var] = d_lambda[0].size() - 1;

      lhs += mult;
      alpha_lhs += problem.d_u1[var] * mult;
    }
    d_constrs[0].push_back(d_cglp.addConstr(lhs <= d_beta[var]));
  }
  d_constrs[0].push_back(d_cglp.addConstr(0 <= d_tau));
  d_constrs[0].push_back(d_cglp.addConstr(0 <= d_kappa));
  d_constrs[0].push_back(d_cglp.addConstr(alpha_lhs >= d_alpha));

  delete[] lambda;
  d_cglp.update();
}