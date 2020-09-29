#include "benders.h"
#include <fstream>

void Benders::reverse_cut(double UB)
{
  d_master.reverse_cut(UB);
  d_agg.reverse_cut(UB);
  d_pslp.reverse_cut(UB);
}