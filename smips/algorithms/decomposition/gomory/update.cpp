#include "gomory.h"

void Gomory::update(vector<double> &rhs, vector<int> &vbasis, vector<int> &cbasis)
{  

  for (size_t con = 0; con != d_ss_leq; ++con)       
  {                               // <= constraints  
    if (cbasis[con] == 0)         // if constraint is basic (non-binding), then
    {
      rhs[con] = GRB_INFINITY;             // relax the constraint by setting rhs = inf
    }
  }

  for (size_t con = d_ss_leq; con != d_ss_leq + d_ss_geq; ++con)       
  {                               // >= constraints  
    if (cbasis[con] == 0)         // if constraint is basic (non-binding), then
      rhs[con] = -GRB_INFINITY;            // relax the constraint by setting rhs = -inf
  }

      // update rhs (implicitly relaxing the appropriate constraints)
  d_model.set(GRB_DoubleAttr_RHS, d_constrs, rhs.data(), rhs.size());

      // relax appropriate variable bounds
  double lb[d_n2];
  double ub[d_n2];
  for (size_t var = 0; var != d_n2; ++var)
  {
    lb[var] = (vbasis[var] == -1) ? d_l2[var] : -GRB_INFINITY;        // lb = -inf if if var not at lower bound
    ub[var] = (vbasis[var] == -2) ? d_u2[var] : GRB_INFINITY;        // ub =  inf if x[var] not at upper bound
  }
    
  d_model.set(GRB_DoubleAttr_LB, d_vars, lb, d_n2);   // setting the lower bounds to their appropriate values  
  d_model.set(GRB_DoubleAttr_UB, d_vars, ub, d_n2);   // setting the lower bounds to their appropriate values 
}
