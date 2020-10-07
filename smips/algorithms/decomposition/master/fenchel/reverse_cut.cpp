#include "fenchel.h"

void Fenchel::reverse_cut(double UB)
{
  if (d_rcut)
  {
    GRBConstr rcut = d_sub.getConstrByName("rcut");
    rcut.set(GRB_DoubleAttr_RHS, UB);
  } else
  {
    GRBLinExpr lhs = d_theta;                        // kappa * theta + beta^T x
    lhs.addTerms(d_problem.d_c.data(), d_xvars.data(), d_xvars.size());
    d_sub.addConstr(lhs, GRB_LESS_EQUAL, UB, "rcut");
    d_rcut = true;
  }

      // removing points (constraints) from mp which violate cx + theta <= UB
  GRBConstr *mp_cons = d_mp.getConstrs();
  for (size_t con = d_points.size() - 1; con != -1; --con)
  {
    Point &point = d_points[con];
    double cx_theta = inner_product(point.d_xvals.begin(), point.d_xvals.end(), d_problem.d_c.begin(), point.d_theta);
    if (cx_theta > UB)
    {
      d_points.erase(d_points.begin() + con);
      d_mp.remove(mp_cons[con]);
    }
  }
  delete[] mp_cons;

  d_sub.update();
  d_mp.update();
}