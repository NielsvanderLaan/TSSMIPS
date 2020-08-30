#include "lagrangian.h"

BendersCut Lagrangian::strong_cut(size_t s, vector<double> &pi)
{
  update(s, pi);
  return BendersCut {solve(), pi, 0.0};
}