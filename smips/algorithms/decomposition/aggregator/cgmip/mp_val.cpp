#include "cgmip.h"

double CGMip::mp_val()
{
  try
  {
    return d_mp.get(GRB_DoubleAttr_ObjVal) + d_alpha.get(GRB_DoubleAttr_X) - d_sub.get(GRB_DoubleAttr_ObjBound);
  } catch (GRBException &e)
  {
    cout << "error in mp_val(), status: " << d_mp.get(GRB_IntAttr_Status) << '\n';
    return GRB_INFINITY;
  }
}