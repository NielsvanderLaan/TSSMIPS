#include "benders.h"
#include <fstream>

void Benders::reverse_cut(double UB)
{
  d_master.reverse_cut(UB);
  d_agg.reverse_cut(UB);
    // TODO delegate to d_pslp
    // d_pslp in turn delegates to d_zk whose elements delegate to their cglps
}