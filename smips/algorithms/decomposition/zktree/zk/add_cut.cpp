#include "zk.h"

#include <cmath>
static bool abs_compare(double a, double b)
{
  return (std::abs(a) < std::abs(b));
}

bool ZK::add_cut(Cut cut, double *x, double theta, double tol, size_t conIdx)
{
  double Tmax = abs(*max_element(cut.Trow.begin(), cut.Trow.end(), abs_compare));
  double Wmax = abs(*max_element(cut.Wrow.begin(), cut.Wrow.end(), abs_compare));
  double abs_max = max(max(Tmax, Wmax), max(abs(cut.r), abs(cut.rhs)));
  if (abs_max < 1e-8)
    return false;

  double scale = 1 / abs_max;

  cut.rhs *= scale;
  cut.r *= scale;
  for_each(cut.Trow.begin(), cut.Trow.end(), [scale](double &val){val *= scale;});
  for_each(cut.Wrow.begin(), cut.Wrow.end(), [scale](double &val){val *= scale;});

    // computing lhs
  double lhs = 0;
  for (size_t var = 0; var != d_n2; ++var)
    lhs += cut.Wrow[var] * d_yvals[var];

    // computing rhs
  double rhs = cut.rhs - cut.r * theta;
  for (size_t var = 0; var != d_n1; ++var)
    rhs -= cut.Trow[var] * x[var];
  
  if (lhs >= rhs - tol)
    return false;
  
    // updating cut coefficients
  d_omega.push_back(cut.rhs);
  d_Wmat.push_back(cut.Wrow);
  d_Tmat.push_back(cut.Trow);
  d_tau.push_back(cut.r);
  d_signs.push_back(1);
  
    // add cut 
  int cind[d_nVars];              // variable indices
  iota(cind, cind + d_n2, 0);  
  GRBaddconstr(d_model, d_n2, cind, cut.Wrow.data(), GRB_EQUAL, rhs, NULL);
  
    // add slack variable
  int vind[1]; vind[0] = conIdx;
  double vval[] = {-1};  
  GRBaddvar(d_model, 1, vind, vval, 0, 0, GRB_INFINITY, GRB_CONTINUOUS, NULL);
  
    //  updating the cglp
  add_cglp_row(cut.Trow.data(), cut.r, cut.Wrow.data(), cut.rhs);

  return true;
}   








