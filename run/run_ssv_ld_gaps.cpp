#include "run.h"

void run_ssv_ld_gaps(Data &rand, GRBEnv &env, GRBenv *c_env)
{
  vector<bool> fs_cont {false, true};
  vector<bool> ss_bin {false, true};
  vector<bool> standard_T {false, true};
  vector<double> size {11};

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
          DeqForm DEF(env, problem);
          DEF.solve(300.0);

          Benders ben(env, c_env, problem, false);
          auto t1 = chrono::high_resolution_clock::now();
          double LP = ben.lpSolve();
          auto t2 = chrono::high_resolution_clock::now();
          double LD = ben.ldSolve();
          auto t3 = chrono::high_resolution_clock::now();
          double SC = ben.ldSolve(false);
          auto t4 = chrono::high_resolution_clock::now();

          cout << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << ' '
               << chrono::duration_cast<chrono::milliseconds>(t3 - t1).count() / 1000.0 << ' '
               << chrono::duration_cast<chrono::milliseconds>(t4 - t1).count() / 1000.0 << '\n';

          cout << DEF.d_objBound << ' ' << DEF.d_objVal << ' ' << LP << ' ' << LD << ' ' << SC << '\n';
        }

      }

    }
  }

}