#include <iostream>
#include <chrono>
#include <vector>

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
      //run_ssv_ld_gaps(rand, env, c_env);
      // create problem
      //Problem problem(10, 0, 0, 5, 5, 5, 100, rand, env, 0, 0, 0, 5);
      //problem.randomInstance();
      //problem.enforce_ccr(50.0);


      Problem problem(rand, env);
      //problem.sslp(15, 45, 5);
      //problem.ssv95(11, 1, 0, 1);
      problem.caroe(100);
      problem.enforce_ccr(1e4);

    /*
      string instance(argv[1]);
      if (instance == "SIZES")
      {
        cout << "SIZES" << argv[2] << endl;
        problem.sizes(stoi(argv[2]));
        problem.enforce_ccr(1e4);
      }
      if (instance == "DCAP")
      {
        cout << "DCAP_" << argv[2] << '_' << argv[3] << '_' << argv[4] << '_' << argv[5] << ' ' << argv[6] << ' ' << argv[7] << '\n';
        problem.dcap(stoi(argv[2]), stoi(argv[3]), stoi(argv[4]), stoi(argv[5]), stoi(argv[6]), stod(argv[7]));
      }
      */



      vector<Type> types = string_to_type(argv, argc);
      for_each(types.begin(), types.end(), [](Type type){cout << name(type) << "s\n";});





      /*
      Tree tree(env, c_env, problem);
      auto t1 = chrono::high_resolution_clock::now();
      vector<double> x_bab = tree.bab(true, false, true, true, false);
      auto t2 = chrono::high_resolution_clock::now();
      for_each(x_bab.begin(), x_bab.end(), [](double val) { cout << val << ' '; });
      cout << "\ncx + Q(x) = " << problem.evaluate(x_bab.data()) << '\n';
      cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';
      */

      /*
      DeqForm DEF(env, problem);
      DEF.d_model.set(GRB_IntParam_OutputFlag, 1);
      DEF.solve(300.0);
      cout << "eta_star = " << DEF.d_objVal << ". LB = " << DEF.d_objBound << '\n';
      */

      {
        auto t1 = chrono::high_resolution_clock::now();
        Benders ben(env, c_env, problem);
        //ben.lpSolve();
        ben.hybrid_solve(types, false, 10000, GRB_INFINITY, 1e-4, 24*3600);
        auto t2 = chrono::high_resolution_clock::now();
        cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';
      }

    }

    GRBfreeenv(c_env);
  } catch (GRBException e)
  {
    cout << e.getErrorCode() << ' ' << e.getMessage() << endl;
  }
}












