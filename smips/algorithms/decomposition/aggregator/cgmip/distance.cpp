#include "cgmip.h"

double CGMip::distance(Point &first, Point &second)
{
  double ret = abs(first.d_eta - second.d_eta);
  ret += abs(first.d_theta - second.d_theta);
  for (size_t idx = 0; idx != first.d_x.size(); ++idx)
    ret += abs(first.d_x[idx] - second.d_x[idx]);

  return ret;
}