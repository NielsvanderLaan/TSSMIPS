#include "zk.h"

bool ZK::optimize()
{
  GRBoptimize(d_model);
  
  int status;
  GRBgetintattr(d_model, "Status", &status);
  if (status == 3)      // model is infeasible
    return false;
  
  GRBgetdblattr(d_model, "ObjVal", &d_objVal);
  GRBgetdblattrarray(d_model, GRB_DBL_ATTR_X, 0, d_n2, d_yvals.data());
  return true;
}