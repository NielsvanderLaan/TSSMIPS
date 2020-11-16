#include "cglp.h"

void Cglp::update_bound(size_t var, double val, bool lower, bool fs)
{
  if (not d_used)
    return;

  vector<int> &mults = fs ? (lower ? d_l1_mults : d_u1_mults) : (lower ? d_l2_mults : d_u2_mults);

  int var_idx = mults[var];
  size_t rhs_con = d_n1 + d_n2 + 1; // index of final constraint (c.t. rhs)
  if (var_idx == -1)            // constraint was of the form x >= 0 or x < 1e10, bound not included.
  {
    size_t var_con = fs ? var : d_n1 + 1 + var;  // constraint index of cglp corresponding to x or y variable
    
    GRBConstr constrs1[2] = {d_constrs1[var_con], d_constrs1[rhs_con]};    // c.t. x[var] and rhs
    GRBConstr constrs2[2] = {d_constrs2[var_con], d_constrs2[rhs_con]};    // idem
    double coeffs[2] = {1, val};                
    if (lower)
    {
      d_lambda1.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, 2, constrs1, coeffs));
      d_lambda2.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, 2, constrs2, coeffs));
    } else
    {
      d_lambda1.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS, 2, constrs1, coeffs));
      d_lambda2.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS, 2, constrs2, coeffs));
    }
    mults[var] = d_nMults;
    ++d_nMults;
  } else
  { 
    d_model.chgCoeff(d_constrs1[rhs_con], d_lambda1[var_idx], val);
    d_model.chgCoeff(d_constrs2[rhs_con], d_lambda2[var_idx], val);
  }
  
  d_model.update();
}