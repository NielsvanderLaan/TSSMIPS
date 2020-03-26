#include "benders.h"

bool Benders::round_of_cuts(Master::Solution sol, double tol)
{
  bool fail = true;
  vector<BendersCut> cuts = d_master.round_of_cuts();

  for (BendersCut cut : cuts)
  {
    cout << cut.d_alpha << '\n';
    for_each(cut.d_beta.begin(), cut.d_beta.end(), [](double val){cout << val << ' ';});
    cout << '\n';
    cout << cut.d_tau << '\n';
  }

  for (BendersCut cut : cuts)
  {
    if (add_cut(cut, sol, tol))
      fail = false;  
  }
  return false;  
}