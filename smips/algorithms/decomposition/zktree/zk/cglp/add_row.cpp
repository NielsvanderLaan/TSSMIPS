#include "cglp.h"

void Cglp::add_row(double *coef_x, double coef_theta, double *coef_y, double rhs)
{
  ++d_nMults;
  
  GRBVar lambda1 = d_model.addVar(0, INFINITY, 0, GRB_CONTINUOUS);
  GRBVar lambda2 = d_model.addVar(0, INFINITY, 0, GRB_CONTINUOUS);
  d_lambda1.push_back(lambda1);
  d_lambda2.push_back(lambda2);
  
  size_t nConstrs = d_n1 + d_n2 + 2;
  GRBVar lambda1_ptr[nConstrs];
  GRBVar lambda2_ptr[nConstrs];
  fill_n(lambda1_ptr, nConstrs, lambda1);
  fill_n(lambda2_ptr, nConstrs, lambda2);
  
  vector<double> vals(nConstrs);
  copy_n(coef_x, d_n1, vals.begin());  
  vals[d_n1] = coef_theta;
  copy_n(coef_y, d_n2, &vals[d_n1 + 1]);
  vals[d_n1 + d_n2 + 1] = rhs;
  
  //cout << "nConstrs = " << nConstrs << '\n';  
  //cout << "numConstrs = " << d_model.get(GRB_IntAttr_NumConstrs) << '\n';
  
  d_model.chgCoeffs(d_constrs1.data(), lambda1_ptr, vals.data(), nConstrs);
  d_model.chgCoeffs(d_constrs2.data(), lambda2_ptr, vals.data(), nConstrs);
  
  
  d_model.update();
}