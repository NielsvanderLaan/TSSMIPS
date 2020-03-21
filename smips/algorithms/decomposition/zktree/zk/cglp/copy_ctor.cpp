#include "cglp.h"

Cglp::Cglp(const Cglp &other):
  d_model(other.d_model),
  d_n1(other.d_n1), 
  d_p1(other.d_p1), 
  d_m1(other.d_m1), 
  d_n2(other.d_n2), 
  d_p2(other.d_p2), 
  d_m2(other.d_m2),
  d_nMults(other.d_nMults),
  d_constrs1(d_n1 + d_n2 + 2),
  d_constrs2(d_n1 + d_n2 + 2),
  d_l1_mults(other.d_l1_mults),
  d_u1_mults(other.d_u1_mults),
  d_l2_mults(other.d_l2_mults),
  d_u2_mults(other.d_u2_mults),
  d_Trow(d_n1),
  d_Wrow(d_n2),
  d_lambda1(d_nMults),
  d_lambda2(d_nMults)
{  
  GRBConstr *constrs = d_model.getConstrs();
  size_t nConstrs = d_n1 + d_n2 + 2;
  for (size_t con = 0; con != nConstrs; ++con)
  {
    d_constrs1[con] = constrs[2 * con];
    d_constrs2[con] = constrs[2 * con + 1];
  }
  delete[] constrs;
  
  GRBVar *vars = d_model.getVars();
  
  for (size_t var = 0; var != d_n1; ++var)
    d_Trow[var] = vars[var];
    
  for (size_t var = 0; var != d_n2; ++var)
    d_Wrow[var] = vars[d_n1 + var];
    
  d_r = vars[d_n1 + d_n2];
  d_h = vars[d_n1 + d_n2 + 1];
  
  size_t begin = d_n1 + d_n2 + 2;
  for (size_t mult = 0; mult != d_nMults; ++mult)
  {
    d_lambda1[mult] = vars[begin + 2*mult];
    d_lambda2[mult] = vars[begin + 2*mult + 1];
  }
  
  delete[] vars;
}