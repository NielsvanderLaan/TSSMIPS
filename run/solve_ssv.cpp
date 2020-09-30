#include "run.h"

void solve_ssv(Data &rand, GRBEnv &env, GRBenv *c_env, vector<Type> types, bool rcuts)
{
  vector<bool> fs_cont {false, true};
  vector<bool> ss_bin {false, true};
  vector<bool> standard_T {false, true};
  vector<double> size {11, 21};

  for (bool fs : {false, true})
  {
    for (bool ss : {false, true})
    {
      for (bool type : {false, true})
      {
        for (double S : {11, 21, 101})
        {
          cout << "SSV_" << S << '_' << fs << '_' << ss << '_' << type << '\n';
          Problem problem(rand, env);
          problem.ssv95(S, fs, ss, type);

          auto t1 = chrono::high_resolution_clock::now();
          Benders ben(env, c_env, problem);
          ben.hybrid_solve(types, false,100000, GRB_INFINITY, 1e-4, 12 * 3600, rcuts);
        }
      }
    }
  }

}