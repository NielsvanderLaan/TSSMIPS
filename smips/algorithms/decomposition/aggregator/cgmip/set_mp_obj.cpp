#include "cgmip.h"

void CGMip::set_mp_obj(double *x, double &theta)
{
  // -S - alpha + beta^Tx + tau * theta = 0
  GRBConstr con = d_mp.getConstrByName("objective");
  d_mp.chgCoeff(con, d_tau, theta);
  vector<GRBConstr> con_vec(d_beta.size(), con);
  d_mp.chgCoeffs(con_vec.data(), d_beta.data(), x,d_beta.size());
  d_mp.update();


  //d_mp.set(GRB_DoubleAttr_Obj, d_beta.data(), x, d_beta.size());
  //d_tau.set(GRB_DoubleAttr_Obj, theta);

}