#include "zktree.h"

void ZkTree::add_child(size_t node_idx)
{
  d_nodes.emplace_back(new ZK(*d_nodes[node_idx]));
  
  d_lb_mult_inds.push_back(d_lb_mult_inds[node_idx]);
  d_ub_mult_inds.push_back(d_ub_mult_inds[node_idx]);
  d_rcut_inds.push_back(d_rcut_inds[node_idx]);

  size_t nMults = d_lambda[node_idx].size();
  double *lb_mults = d_cglp.get(GRB_DoubleAttr_LB, d_lambda[node_idx].data(), nMults);
  double *ub_mults = d_cglp.get(GRB_DoubleAttr_UB, d_lambda[node_idx].data(), nMults);

  string names[nMults];
  string base = "lambda_" + to_string(d_nodes.size() - 1) + "_";
  for (size_t mult = 0; mult != nMults; ++mult)
    names[mult] = base + to_string(mult);
  GRBVar *mults = d_cglp.addVars(lb_mults, ub_mults, NULL, NULL, names, nMults);
  d_lambda.push_back(vector<GRBVar>(mults, mults + nMults));
  delete[] lb_mults, ub_mults;

  vector<GRBConstr> constrs;

  double coeff[nMults];
  for (size_t idx = 0; idx != d_beta.size(); ++idx)
  {
    GRBLinExpr lhs;
    for (size_t mult = 0; mult != nMults; ++mult)
      lhs += d_cglp.getCoeff(d_constrs[node_idx][idx], d_lambda[node_idx][mult]) * mults[mult];
    constrs.push_back(d_cglp.addConstr(lhs <= d_beta[idx]));
  }

  size_t idx = d_beta.size();
  for (size_t mult = 0; mult != nMults; ++mult)
    coeff[mult] = d_cglp.getCoeff(d_constrs[node_idx][idx], d_lambda[node_idx][mult]);
  GRBLinExpr tau_lhs;
  tau_lhs.addTerms(coeff, mults, nMults);
  constrs.push_back(d_cglp.addConstr(tau_lhs <= d_tau));

  ++idx;
  for (size_t mult = 0; mult != nMults; ++mult)
    coeff[mult] = d_cglp.getCoeff(d_constrs[node_idx][idx], d_lambda[node_idx][mult]);
  GRBLinExpr kappa_lhs;
  kappa_lhs.addTerms(coeff, mults, nMults);
  constrs.push_back(d_cglp.addConstr(kappa_lhs <= d_kappa));

  ++idx;
  for (size_t mult = 0; mult != nMults; ++mult)
    coeff[mult] = d_cglp.getCoeff(d_constrs[node_idx][idx], d_lambda[node_idx][mult]);
  GRBLinExpr alpha_lhs;
  alpha_lhs.addTerms(coeff, mults, nMults);
  constrs.push_back(d_cglp.addConstr(alpha_lhs >= d_alpha));

  d_constrs.push_back(constrs);

  d_cglp.update();
  delete[] mults;
}