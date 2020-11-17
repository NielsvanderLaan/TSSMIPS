#include "zktree.h"

void ZkTree::add_row_to_cglp(const double *coeff_x, double coeff_theta, double coeff_eta, double rhs, size_t node_idx)
{
  string name = "lambda_" + to_string(node_idx) + "_" + to_string(d_lambda[node_idx].size());
  GRBVar lambda = d_cglp.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
  d_lambda[node_idx].push_back(lambda);

  vector<GRBVar> lambda_vec(d_beta.size(), lambda);
  d_cglp.chgCoeffs(d_constrs[node_idx].data(), lambda_vec.data(), coeff_x, d_beta.size());

  size_t con_idx = d_beta.size();
  d_cglp.chgCoeff(d_constrs[node_idx][con_idx], lambda, coeff_theta);
  d_cglp.chgCoeff(d_constrs[node_idx][con_idx + 1], lambda, coeff_eta);
  d_cglp.chgCoeff(d_constrs[node_idx][con_idx + 2], lambda, rhs - (coeff_eta + coeff_theta) * d_L);

  d_cglp.update();
}
