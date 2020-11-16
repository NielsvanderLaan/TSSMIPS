#include "cglp.h"

void Cglp::create_disjunction(size_t var_idx, double val)
{
  if (not d_used)
    return;

      // pointer to constraints of first disjunction corresponding to the y variables
  GRBConstr *constrs1 = &d_constrs1[d_n1 + 1];  
  GRBConstr *constrs2 = &d_constrs2[d_n1 + 1];  // idem for second disjunction
  
  vector<GRBVar> vars1(d_n2, d_lambda1[0]);  
  vector<GRBVar> vars2(d_n2, d_lambda2[0]); 
  
  vector<double> vals(d_n2);
  vals[var_idx] = 1.0;

  d_model.chgCoeffs(constrs1, vars1.data(), vals.data(), d_n2);
  d_model.chgCoeffs(constrs2, vars2.data(), vals.data(), d_n2);
  
  d_model.chgCoeff(d_constrs1[d_n1 + d_n2 + 1], d_lambda1[0], val);
  d_model.chgCoeff(d_constrs2[d_n1 + d_n2 + 1], d_lambda2[0], val + 1);
  
  d_model.update();
}