#include "run.h"

void solve_dcap(bool lp_cuts, bool sb_cuts, bool zk_cuts, bool strong_cuts, Data &rand, GRBEnv &env, GRBenv *c_env)
{
  vector<size_t> nsamples {200, 300, 500};
  vector<vector<size_t>> inputs {{2, 3, 3},
                                 {2, 4, 3},
                                 {3, 3, 2},
                                 {3, 4, 2}};
  for (vector<size_t> size : inputs)
  {
    for (size_t S : nsamples )
    {
      cout << "INSTANCE: " << size[0] << ' ' << size[1] << ' ' << size[2] << ' ' << S << '\n';

      Problem problem(rand, env);
      problem.dcap(size[0], size[1], size[2], S);
      Tree tree(env, c_env, problem);
      auto t1 = chrono::high_resolution_clock::now();
      vector<double> x_bab = tree.bab(lp_cuts, sb_cuts, zk_cuts, strong_cuts, false, true);
      auto t2 = chrono::high_resolution_clock::now();
      cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << endl;


    }

  }

}