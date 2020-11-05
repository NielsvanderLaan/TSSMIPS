#include "zk.h"
#include <cmath>

Cut ZK::generate_gmi_cut(Master &master, size_t row, double yval, double *x, bool zk, bool check)
{
  GRBmodel *model = master.d_cmodel;
  vector<double> &kappa = master.d_kappa;
  vector<vector<double>> &beta = master.d_beta;
  vector<double> &gamma = master.d_gamma;

  int nVarsMaster;                                      // number of variables in master problem (including slacks)
  GRBgetintattr(model, "NumVars", &nVarsMaster);

      // computing tableau row
  double tab_row_x[nVarsMaster];
  compute_tab_row_x(tab_row_x, nVarsMaster, row, model, zk); // tableau row for (theta, x) (in that order)

  double tab_row_y[d_nVars];
  compute_tab_row_y(tab_row_y, row); // tableau row for y variables 
  
  double coef_x[nVarsMaster - 1]; double coef_y[d_nVars]; double coef_theta = -1; double coef_rhs = 1;// cut coefficients
  double a0 = zk ? yval : yval + inner_product(x, x + d_n1, tab_row_x + 1, 0.0);

  double mastervars[nVarsMaster];
  GRBgetdblattrarray(master.d_cmodel, "X", 0, nVarsMaster, mastervars);
  double yvals[d_nVars];
  GRBgetdblattrarray(d_model, "X", 0, d_nVars, yvals);
  a0 = inner_product(tab_row_x, tab_row_x + nVarsMaster, mastervars, 0.0) + inner_product(tab_row_y, tab_row_y + d_nVars, yvals, 0.0);


  bool proper = gmi_cut(tab_row_x, tab_row_y, a0, coef_x, coef_y, coef_theta, nVarsMaster);
  if (not proper)
    return Cut{ vector<double>(d_n1, 0.0), 0.0, vector<double>(d_n2, 0.0), 0.0};

  transform_cut(coef_x, coef_y, coef_theta, coef_rhs, kappa, beta, gamma, nVarsMaster - d_n1 - 1);
  vector<double> Trow(coef_x, coef_x + d_n1);
  vector<double> Wrow(coef_y, coef_y + d_n2);

  return Cut { Trow, coef_theta, Wrow, coef_rhs };      
}