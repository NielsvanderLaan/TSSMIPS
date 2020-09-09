#include "benders.h"

size_t Benders::round_of_cuts(Master::Solution sol, double tol)
{
  size_t count = 0;

  vector<BendersCut> cuts = d_master.round_of_cuts();

  for (BendersCut cut : cuts)
  {
    if (not add_cut(cut, sol, tol))
      ++count;
  }

  return count;
}