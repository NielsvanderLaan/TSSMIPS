#include "lagrangian.h"

void Lagrangian::reverse_cut(double UB)
{
  if (d_rcut)
  {
    GRBConstr rcut = d_model.getConstrByName("rcut");
    rcut.set(GRB_DoubleAttr_RHS, UB);
  } else
  {
    GRBLinExpr lhs = d_theta;
    lhs.addTerms(d_problem.d_c.data(), d_z_vars.data(), d_z_vars.size());
    d_model.addConstr(lhs <= UB, "rcut");
    d_rcut = true;
  }

  d_model.update();
}