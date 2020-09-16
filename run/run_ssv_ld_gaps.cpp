#include "run.h"

void run_ssv_ld_gaps(Data &rand, GRBEnv &env, GRBenv *c_env)
{
  vector<bool> fs_cont {false, true};
  vector<bool> ss_bin {false, true};
  vector<bool> standard_T {false, true};
  vector<double> size {11, 21};

  for (bool fs : fs_cont)
  {
    for (bool ss : ss_bin)
    {
      for (bool type : standard_T)
      {
        for (double S : size)
        {
          cout << "INSTANCE: " << fs << ' ' << ss << ' ' << type << ' ' << S << '\n';
          Problem problem(rand, env);
          problem.ssv95(S, fs, ss, type);

          auto t1 = chrono::high_resolution_clock::now();
          Benders ben(env, c_env, problem);
          ben.hybrid_solve(true, false, false, true, false, false, 10000, GRB_INFINITY);
          auto t2 = chrono::high_resolution_clock::now();
          cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';
        }

      }

    }
  }

}