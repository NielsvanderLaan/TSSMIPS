#include "cglp.h"

Cut Cglp::generate_cut(double *x, double theta, double *y, size_t var_idx, double val)
{
  create_disjunction(var_idx, val);
  set_obj(x, theta - d_L, y);    // cglp is in terms of theta' = theta - L

  d_model.set(GRB_IntParam_NumericFocus, 0);
  d_model.optimize();

  int status = d_model.get(GRB_IntAttr_Status);
  if (status != 2)
  {
    d_model.set(GRB_IntParam_NumericFocus, 3);
    d_model.reset();
    d_model.optimize();
    status = d_model.get(GRB_IntAttr_Status);
    if (status != 2)
    {
      d_model.write("cg.lp");
      exit(status);
    }
  }

  double violation = d_model.get(GRB_DoubleAttr_ConstrVio);
  double resid = d_model.get(GRB_DoubleAttr_ConstrResidual);
  if (violation + resid > 1e-6)
  {
    cout << "violation: " << violation << ". residual: " << resid <<'\n';
  }




  double r = d_r.get(GRB_DoubleAttr_X);
  double rhs = d_h.get(GRB_DoubleAttr_X) + r * d_L;   // cglp is in terms of theta' = theta - L
  
  double *Trow_values = d_model.get(GRB_DoubleAttr_X, d_Trow.data(), d_n1);
  double *Wrow_values = d_model.get(GRB_DoubleAttr_X, d_Wrow.data(), d_n2);
  
  vector<double> Trow(Trow_values, Trow_values + d_n1);
  vector<double> Wrow(Wrow_values, Wrow_values + d_n2);
  
  delete[] Trow_values;
  delete[] Wrow_values;

          // strengthening the cut
  double lambda1_0 = d_lambda1[0].get(GRB_DoubleAttr_X); 
  double lambda2_0 = d_lambda2[0].get(GRB_DoubleAttr_X); 
  double diff = lambda1_0 - lambda2_0;
  
  if (diff > -1e-8)
    return Cut{ Trow, r, Wrow, rhs };  
    
  for (size_t var = 0; var != d_p1; ++var)
  {
    double coeff = Trow[var];
          // in gurobi slack := right - left. We need to know left = lambda_aj, where right = coeff    
    double lambda1_aj = coeff - d_constrs1[var].get(GRB_DoubleAttr_Slack);
    double lambda2_aj = coeff - d_constrs2[var].get(GRB_DoubleAttr_Slack);
  
    double m = (lambda2_aj - lambda1_aj) / diff;
  
    Trow[var] = min(lambda1_aj + lambda1_0 * floor(m), lambda2_aj + lambda2_0 * ceil(m));
  }
  
  for (size_t var = 0; var != d_p2; ++var)
  {
    if (var == var_idx)
      continue;
      
    double coeff = Wrow[var];
    size_t con = d_n1 + 1 + var;   // first d_n1 + 1 constraints correspond to x and theta, we need the y variables
    double lambda1_aj = coeff - d_constrs1[con].get(GRB_DoubleAttr_Slack);
    double lambda2_aj = coeff - d_constrs2[con].get(GRB_DoubleAttr_Slack);
  
    double m = (lambda2_aj - lambda1_aj) / diff;
  
    Wrow[var] = min(lambda1_aj + lambda1_0 * floor(m), lambda2_aj + lambda2_0 * ceil(m));
  }

  return Cut { Trow, r, Wrow, rhs };
}










