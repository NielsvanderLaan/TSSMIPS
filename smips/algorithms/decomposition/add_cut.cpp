#include "benders.h"
#include <fstream>

bool Benders::add_cut(BendersCut &cut, Master::Solution sol, double tol)
{
  if (cut.d_tau > 0)
    cut.scale();

  bool improper_cut = d_master.add_cut(cut, sol, tol);

  if (not improper_cut)
    d_agg.add_rows(cut);

  return improper_cut;
}