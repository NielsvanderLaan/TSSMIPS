#include "benders.h"

void Benders::update(double UB)
{
  d_UB = UB;
  reverse_cut(UB);
}