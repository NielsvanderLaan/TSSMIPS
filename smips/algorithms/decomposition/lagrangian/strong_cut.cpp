#include "lagrangian.h"

BendersCut Lagrangian::strong_cut(vector<double> &pi)
{
  update(pi);
  return BendersCut {solve(), pi, 0.0};
}