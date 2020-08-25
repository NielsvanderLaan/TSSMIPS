#include "run.h"

void solve_sizes(Data &rand, GRBEnv &env, GRBenv *c_env)
{
  vector<size_t> nsamples {3, 5, 10};

  for (size_t S : nsamples )
  {
    cout << "INSTANCE: " << S << '\n';

    Problem problem(rand, env);
    problem.sizes(S);
    problem.enforce_ccr(1e4);
    Tree tree(env, c_env, problem);
    auto t1 = chrono::high_resolution_clock::now();
    vector<double> x_bab = tree.bab( false, 1e-2);
    auto t2 = chrono::high_resolution_clock::now();
    cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << endl;

  }


}
