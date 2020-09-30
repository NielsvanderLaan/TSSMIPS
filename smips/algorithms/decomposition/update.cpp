#include "benders.h"

void Benders::update(double UB, bool rcuts)
{
  d_UB = UB;
  if (rcuts) reverse_cut(UB);
}