#include "cglp.h"

void Cglp::add_row(double *coef_x, double coef_theta, double *coef_y, double rhs, bool geq)
{
  if (not d_used)
    return;

  ++d_nMults;

  GRBVar lambda1 = d_model.addVar(geq ? 0 : -GRB_INFINITY, geq ? INFINITY : 0, 0, GRB_CONTINUOUS);
  GRBVar lambda2 = d_model.addVar(geq ? 0 : -GRB_INFINITY, geq ? INFINITY : 0, 0, GRB_CONTINUOUS);
  d_lambda1.push_back(lambda1);
  d_lambda2.push_back(lambda2);
  
  size_t nConstrs = d_n1 + d_n2 + 2;
  vector<GRBVar> lambda1_vec(nConstrs, lambda1);
  vector<GRBVar> lambda2_vec(nConstrs, lambda2);

  vector<double> vals(nConstrs);
  copy_n(coef_x, d_n1, vals.begin());  
  vals[d_n1] = coef_theta;
  copy_n(coef_y, d_n2, &vals[d_n1 + 1]);
  vals[d_n1 + d_n2 + 1] = rhs - coef_theta*d_L;         // cglp is in terms of theta' = theta - L



  d_model.chgCoeffs(d_constrs1.data(), lambda1_vec.data(), vals.data(), nConstrs);
  d_model.chgCoeffs(d_constrs2.data(), lambda2_vec.data(), vals.data(), nConstrs);

  d_model.update();
}