#include <iostream>
#include <chrono>
#include <vector>

#include "gurobi_c++.h"
#include "gurobi_c.h"

#include "smips/problem_data/problem.h"
#include "smips/algorithms/deqform/deqform.h"
#include "smips/algorithms/decomposition/benders.h"
#include "smips/algorithms/trees/tree.h"
#include "smips/algorithms/decomposition/zktree/zktree.h"

#include "run/run.h"

using namespace std;

int main(int argc, char *argv[])
{
  try
  {
    Data rand(14785);  // test

    GRBEnv env;
    env.set(GRB_IntParam_OutputFlag, 0);
    env.set(GRB_IntParam_Threads, 1);

    GRBenv *c_env;
    GRBloadenv(&c_env, nullptr);
    GRBsetintparam(c_env, "OutputFlag", 0);
    GRBsetintparam(c_env, "Threads", 1);

    {
      //solve_ri(rand, env, c_env);
      // create problem
      //Problem problem(10, 0, 0, 5, 5, 5, 100, rand, env, 0, 0, 0, 5);
      //problem.randomInstance();
      //problem.enforce_ccr(50.0);
      Problem problem(rand, env);
      bool sizes = stoi(argv[1]);
      if (sizes)
      {
        cout << "SIZES" << argv[2] << endl;
        problem.sizes(stoi(argv[2]));
        problem.enforce_ccr(1e4);
      } else
      {
        cout << "DCAP_" << argv[2] << '_' << argv[3] << '_' << argv[4] << '_' << argv[5] << '\n';
        problem.dcap(stoi(argv[2]),stoi(argv[3]), stoi(argv[4]),stoi(argv[5]));
      }



      //problem.ssv95(11, 1,1, 0);

      //problem.sslp(15, 45, 5);
      //problem.dcap(2,3,3,200);
      //problem.enforce_ccr(1e4);
      //size_t S = 100;
      //problem.caroe(S);
      //problem.enforce_ccr(1e4);



    /*
      Tree tree(env, c_env, problem);
      auto t1 = chrono::high_resolution_clock::now();
      vector<double> x_bab = tree.bab(true, true, false, false, false);
      auto t2 = chrono::high_resolution_clock::now();
      for_each(x_bab.begin(), x_bab.end(), [](double val) { cout << val << ' '; });
      cout << "\ncx + Q(x) = " << problem.evaluate(x_bab.data()) << '\n';
      cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';
    */

    /*
      DeqForm DEF(env, problem);
      DEF.solve(300.0);
      cout << "eta_star = " << DEF.d_objVal << '\n';

      Problem ld(rand, env);
      ld.caroe_LD(S);
      DeqForm DEF2(env, ld);
      DEF2.solve(300.0);
      cout << "LD = " << DEF2.d_objVal<< '\n';
    */
      {
        auto t1 = chrono::high_resolution_clock::now();
        Benders ben(env, c_env, problem);
        ben.lpSolve();
        ben.hybrid_solve(true, true, false, true, true, false, 10000, GRB_INFINITY);
        auto t2 = chrono::high_resolution_clock::now();
        cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';
      }
      {
        auto t1 = chrono::high_resolution_clock::now();
        Benders ben(env, c_env, problem);
        ben.lpSolve();
        ben.hybrid_solve(true, false, true, true, false, false, 10000, GRB_INFINITY);
        auto t2 = chrono::high_resolution_clock::now();
        cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';
      }








      //cout << "x = " << *ben.d_incumbent << '\n';


      /*
      double LD = ben.ldSolve(true);
      double SC = ben.ldSolve(false);
      cout << "LD = " << LD << ". SC = " << SC << '\n';


      cout << "L-shaped LB: " << lpLB << '\n';
      x = ben.d_xvals;
      cout << "x: ";
      for (size_t var = 0; var != problem.d_n1; ++var)
        cout << x[var] << ' ';
      cout << '\n';
      cout << "cx + Q(x) = " << problem.evaluate(x) << '\n';
      */

    }

    GRBfreeenv(c_env);
  } catch (GRBException e)
  {
    cout << e.getErrorCode() << ' ' << e.getMessage() << endl;
  }
}












