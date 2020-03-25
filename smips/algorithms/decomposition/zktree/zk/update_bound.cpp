#include "zk.h"

void ZK::update_bound(int var, double val, bool lower, bool fs)
{
  d_cglp.update_bound(var, val, lower, fs);
  
  if (fs)      // if first-stage bound, then we are done
    return;
  
  int con = lower ? d_lb_inds[var] : d_ub_inds[var];    // constraint idx of bound (-1 if it does not exist)

  if (con != -1)        // if bound was added previously, only update rhs
  {   
    d_omega[con] = val; 
    GRBsetdblattrelement(d_model, "RHS", con, val);  
    GRBupdatemodel(d_model);
    return;
  }
        
        // bound was not added previously, we add it now  
  d_Wmat.push_back(vector<double>(d_n2));  // update constraint data
  d_Wmat[d_Wmat.size() - 1][var] = 1;      // idem
  d_Tmat.push_back(vector<double>(d_n1));  // idem
  d_tau.push_back(0);                      // idem
  d_omega.push_back(val);                  // idem 
  
  if (lower)                            
  {
    d_signs.push_back(1);                   // >= constraint  
    d_lb_inds[var] = d_nConstrs;            // index of constraint corresponding to bound (for future valls to update_bound()) 
  } else
  {
    d_signs.push_back(-1);
    d_ub_inds[var] = d_nConstrs;
  }

  GRBaddvar(d_model, 0, NULL, NULL, 0.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS, NULL); // slack variable
  int inds[2] = {var, d_nVars};                                            // variable indices
  double vals[2] = {1.0, 1.0};                                             // coefficients
  if (lower)            // x >= l (so slack features with -1.0)
    vals[1] = -1.0;
  GRBaddconstr(d_model, 2, inds, vals, GRB_EQUAL, val, NULL);              // add constraint
  
  ++d_nConstrs;                            // added one constraint
  ++d_nVars;                               // using one slack variable 
  GRBupdatemodel(d_model);
}