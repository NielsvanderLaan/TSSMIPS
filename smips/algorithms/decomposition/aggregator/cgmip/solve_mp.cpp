#include "cgmip.h"

void CGMip::solve_mp()
{
  d_mp.optimize();
  int status = d_mp.get(GRB_IntAttr_Status);
  if (status != 2)
  {
    cout << "mp status: " << status << '\n';
    exit(1);

  }
}