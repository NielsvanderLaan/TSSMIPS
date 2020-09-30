#include "lagrangian.h"

vector<double> Lagrangian::z_vals()
{
  if (d_model.get(GRB_IntAttr_Status) != 2)
  {
    cout << "error in Lagrangian::z_vals(): NOT OPTIMAL\n";
    exit(1);
  }

  double *z_vals = d_model.get(GRB_DoubleAttr_X, d_z_vars.data(), d_n1);
  vector<double> ret(z_vals, z_vals + d_n1);
  delete[] z_vals;

  return ret;
}