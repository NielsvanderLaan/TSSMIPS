#include "cglp.h"

void Cglp::set_obj(double *x, double theta, double *y)
{
  d_model.set(GRB_DoubleAttr_Obj, d_Trow.data(), x, d_n1);
  d_model.set(GRB_DoubleAttr_Obj, d_Wrow.data(), y, d_n2);
  
  if (theta < 1e10)                    
  {
    d_r.set(GRB_DoubleAttr_UB, 1);        
    d_r.set(GRB_DoubleAttr_Obj, theta);
  } else
    d_r.set(GRB_DoubleAttr_UB, 0);

  
  d_model.update();
}