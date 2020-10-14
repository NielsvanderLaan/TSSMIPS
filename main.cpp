#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>

#include "gurobi_c++.h"
#include "gurobi_c.h"

#include "smips/problem_data/problem.h"
#include "smips/algorithms/deqform/deqform.h"
#include "smips/algorithms/decomposition/benders.h"
#include "smips/algorithms/trees/tree.h"

#include "run/run.h"

using namespace std;

int main(int argc, char *argv[])
{
  try
  {
    Data rand(14785);

    GRBEnv env;
    env.set(GRB_IntParam_OutputFlag, 0);
    env.set(GRB_IntParam_Threads, 1);

    GRBenv *c_env;
    GRBloadenv(&c_env, nullptr);
    GRBsetintparam(c_env, "OutputFlag", 0);
    GRBsetintparam(c_env, "Threads", 1);
    {
      Problem problem(rand, env);
      instance(problem, argc, argv);
      vector<Type> types = string_to_type(argc, argv);
      bool rcuts = use_rcuts(argc, argv);
      bool fenchel = use_fenchel(argc, argv);
      int max_rounds = get_max_rounds(argc, argv);
      double time_limit = get_time_limit(argc, argv);
      details(types, max_rounds, rcuts, fenchel, time_limit);

      if (solve_DEF(argc, argv))
      {
        DeqForm DEF(env, problem);
        DEF.d_model.set(GRB_IntParam_OutputFlag, 1);
        DEF.solve(time_limit);
      }

      if (solve_tree(argc, argv))
      {
        Tree tree(env, c_env, problem);
        vector<double> x_bab = tree.bab(types, rcuts, fenchel, max_rounds, 1e-4,time_limit);
      }

      if (solve_root(argc, argv))
      {
        Benders ben(env, c_env, problem, true);
        ben.lpSolve();
        ben.hybrid_solve(types, false, max_rounds, GRB_INFINITY, 1e-4, time_limit, rcuts, fenchel);
      }
    }

    GRBfreeenv(c_env);
  } catch (GRBException e)
  {
    cout << e.getErrorCode() << ' ' << e.getMessage() << endl;
  }
}












