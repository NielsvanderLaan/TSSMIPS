#include "benders.h"

void Benders::report(double time, double LB, double eval_time, size_t evaluations, size_t gmi_cuts, size_t round,
                     double gmi_time, size_t fenchel_cuts, double fenchel_time, vector<Type> &types,
                     vector<double> &times, vector<size_t> nCuts)
{
  cout << "====================================================\n";
  cout << "LB: " << LB << " UB: " << d_UB << '\n';
  cout << "computation time: " << time << '\n';
  cout << "evaluations: " << evaluations << " (" << avg(eval_time, evaluations) << "s)\n";
  cout << "GMI cuts: " << gmi_cuts << " (" << round << " rounds, " << avg(gmi_time, round) << "s)\n";
  cout << "Fenchel cuts: " << fenchel_cuts << " (" << avg(fenchel_time, fenchel_cuts) << "s)\n";
  for (size_t idx = 0; idx != types.size(); ++idx)
    cout << name(types[idx]) << "s: " << nCuts[idx] << " (" << avg(times[idx], nCuts[idx]) << "s)\n";
  cout << "====================================================\n";
}