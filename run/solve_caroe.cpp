#include "run.h"

void solve_caroe(Data &rand, GRBEnv &env, GRBenv *c_env)
{
  vector<Type> types {LR_LAP, LR, SC_RG, SC_LAP, SC_ZK, SC_BAB};

  Problem problem(rand, env);
  problem.caroe(100);
  problem.enforce_ccr(1e4);

  for (Type type : types)
  {
    for (bool rcuts : vector<bool> {false, true})
    {
      cout << name(type) << "s\n";
      if (rcuts)
        cout << "reverse cuts\n";
      Benders ben(env, c_env, problem);
      ben.hybrid_solve(vector<Type> {type}, false, 0, GRB_INFINITY, 1e-4, 1e100, rcuts);
      cout << '\n';
    }
  }
}
