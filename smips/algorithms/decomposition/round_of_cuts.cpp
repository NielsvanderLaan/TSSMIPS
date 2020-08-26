#include "benders.h"



bool Benders::round_of_cuts(Master::Solution sol, double tol)
{
  cout << "Benders::round_of_cuts()" << endl;

  bool fail = true;

  vector<BendersCut> cuts = d_master.round_of_cuts();

  for (BendersCut cut : cuts)
  {
    if (not add_cut(cut, sol, tol))
      fail = false;
  }

  return fail;
}