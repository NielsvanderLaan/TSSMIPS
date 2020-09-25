#include "pslp.h"

void Pslp::reverse_cut(double UB)
{
  for (size_t s = 0; s != d_S; ++s)
    d_zk[s].reverse_cut(UB);
}