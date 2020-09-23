#include "cgmip.h"

void CGMip::sub_optimal()
{
  int status = d_sub.get(GRB_IntAttr_Status);
  double violation = d_sub.get(GRB_DoubleAttr_ConstrVio);
  if (violation > 1e-4)
    cout << "sub violation = " << violation << '\n';

  if (status == 9 or status == 13)
  {
    cout << "sub gap = " << d_sub.get(GRB_DoubleAttr_ObjVal) - d_sub.get(GRB_DoubleAttr_ObjBound) << '\n';
    return;
  }
  if (status != 2)
  {
    cout << "sub status: " << status << '\n';
    d_sub.write("sub.lp");
    exit(1);
  }
}