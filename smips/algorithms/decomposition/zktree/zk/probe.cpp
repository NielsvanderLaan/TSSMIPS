#include "zk.h"

double ZK::probe(size_t var_idx, double val, bool lower)
{
  int con_idx = lower ? d_lb_inds[var_idx] : d_ub_inds[var_idx];  // constraint idx of upper/lower bound
  double obj;
  
  if (con_idx == -1)
  {  
    string bound_type = lower ? "LB" : "UB";  
    const char *attr = bound_type.c_str();
      
    double old_val;
    GRBgetdblattrelement(d_model, attr, var_idx, &old_val);  // store old bound
    GRBsetdblattrelement(d_model, attr, var_idx, val);       // change bound (temporarily)
    GRBoptimize(d_model); 
    
    int status;                                              // retrieve status    
    GRBgetintattr(d_model, "Status", &status);
    if (status == 3)                                         // check if infeasible 
      obj = GRB_INFINITY;                                    // infinite increase in objval
    else
      GRBgetdblattr(d_model, "ObjVal", &obj);                // retrieve objective value
      
    GRBsetdblattrelement(d_model, attr, var_idx, old_val);   // restore bound
  } else
  {
    GRBsetdblattrelement(d_model, "RHS", con_idx, val);      // change bound (temporarily)
    GRBoptimize(d_model);                                    // optimize the model
    
    int status;                                              // retrieve status    
    GRBgetintattr(d_model, "Status", &status);
    if (status == 3)                                         // check if infeasible 
      obj = GRB_INFINITY;                                    // infinite increase in objval
    else
      GRBgetdblattr(d_model, "ObjVal", &obj);       // retrieve objective value
                                                             
    GRBsetdblattrelement(d_model, "RHS", con_idx, d_omega[con_idx]);    // restore bound
  }
  
  GRBupdatemodel(d_model);        // make changes (restore model)
  
  return obj - d_objVal;
}