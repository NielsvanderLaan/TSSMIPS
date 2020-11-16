#include "zktree.h"

void ZkTree::cglp_bounds(double M, bool set)
{
  vector<double> lb(d_beta.size(), -M);
  vector<double> ub(d_beta.size(), M);
  d_cglp.set(GRB_DoubleAttr_LB, d_beta.data(), lb.data(), d_beta.size());
  d_cglp.set(GRB_DoubleAttr_UB, d_beta.data(), ub.data(), d_beta.size());
  d_tau.set(GRB_DoubleAttr_UB, M);

  string l = "left";
  string r = "right";

  if (set)
  {
    GRBLinExpr alpha_prime = d_alpha + d_tau * d_L + d_L;
    d_cglp.addConstr(alpha_prime >= -M, l);
    d_cglp.addConstr(alpha_prime <= M, r);
  } else
  {
    d_cglp.remove(d_cglp.getConstrByName(l));
    d_cglp.remove(d_cglp.getConstrByName(r));
  }

  d_cglp.update();
}