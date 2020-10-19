#include "cgmip.h"

bool CGMip::mp_optimal()
{
  int status =  d_mp.get(GRB_IntAttr_Status);

  if (status == 9)
    cout << "CGMP timed out" << endl;

  return status == 2;
}