#include "fenchel.h"

void Fenchel::set_mp_obj(vector<double> &x, double &theta)
{
  d_kappa.set(GRB_DoubleAttr_Obj, theta);
  d_mp.set(GRB_DoubleAttr_Obj, d_beta.data(), x.data(), d_beta.size());

  d_mp.update();
}