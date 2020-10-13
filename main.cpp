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
      string instance(argv[1]);
      vector<Type> types = string_to_type(argv, argc);
      for_each(types.begin(), types.end(), [](Type type){cout << name(type) << "s\n";});
      bool rcuts = true;
      bool fenchel = true;
      size_t max_rounds = 0;
      for (size_t idx = 1; idx != argc; ++idx)
      {
        string arg(argv[idx]);
        if (arg == "OFF")
          rcuts = false;
        if (arg == "GOMORY")
          fenchel = false;
        string mr("MAXROUNDS=");
        if (arg.find(mr) == 0)
          max_rounds = stoi(arg.substr(mr.size(), arg.size() - mr.size()));
      }
      cout << "rcuts: " << (rcuts ? "yes\n" : "no\n");
      cout << (fenchel ? "Fenchel" : "Gomory") << " cuts\n";
      cout << "max rounds = " << max_rounds << '\n';

      Problem problem(rand, env);

      if (instance == "SIZES")
      {
        cout << "SIZES" << argv[2] << endl;
        problem.sizes(stoi(argv[2]));
        problem.enforce_ccr(1e4);
      }
      if (instance == "DCAP")
      {
        cout << "DCAP_" << argv[2] << '_' << argv[3] << '_' << argv[4] << '_' << argv[5] << ' ' << argv[6] << '\n';
        problem.dcap(stoi(argv[2]), stoi(argv[3]), stoi(argv[4]), stoi(argv[5]), stoi(argv[6]));
      }
      if (instance == "SSV")
      {
        cout << "SSV_" << argv[2] << '_' << argv[3] << '_' << argv[4] << '_' << argv[5] << '\n';
        problem.ssv95(stoi(argv[2]), stoi(argv[3]), stoi(argv[3]), stoi(argv[4]));
      }
      if (instance == "CAROE")
        solve_caroe(rand, env, c_env);

      /*
      DeqForm DEF(env, problem);
      DEF.d_model.set(GRB_IntParam_OutputFlag, 1);
      DEF.solve(7200.0);
      cout << "eta_star = " << DEF.d_objVal << ". LB = " << DEF.d_objBound << '\n';
      for_each(DEF.d_xVals, DEF.d_xVals + problem.d_n1, [](double val){cout << val << ' ';});
      cout << '\n';
      */

      Tree tree(env, c_env, problem);
      vector<double> x_bab = tree.bab(types, rcuts, fenchel, max_rounds, 1e-4,12*3600);

      /*
      Benders ben(env, c_env, problem, true);
      ben.lpSolve();
      ben.hybrid_solve(types, false, 100000, GRB_INFINITY, 1e-4, 12 * 3600, rcuts, fenchel);
      */
    }

    GRBfreeenv(c_env);
  } catch (GRBException e)
  {
    cout << e.getErrorCode() << ' ' << e.getMessage() << endl;
  }
}












