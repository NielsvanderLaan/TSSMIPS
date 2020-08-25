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
    Data rand(51615);  // test

    GRBEnv env;
    env.set(GRB_IntParam_OutputFlag, 0);
    env.set(GRB_IntParam_Threads, 1);

    GRBenv *c_env;
    GRBloadenv(&c_env, nullptr);
    GRBsetintparam(c_env, "OutputFlag", 0);
    GRBsetintparam(c_env, "Threads", 1);

    {
      //solve_sizes(rand, env, c_env);
      //solve_dcap(rand, env, c_env);
      // create problem
      //Problem problem(10, 0, 0, 5, 5, 5, 100, rand, env, 0, 0, 0, 5);
      //problem.randomInstance();
      //problem.enforce_ccr(50.0);

      Problem problem(rand, env);
      //problem.ssv95(11, 0,1, 1);
      problem.sizes(3);

      //problem.sslp(15, 45, 5);
      //problem.dcap(3,4,2,500);
      //problem.enforce_ccr(100);


      Tree tree(env, c_env, problem);
      auto t1 = chrono::high_resolution_clock::now();
      vector<double> x_bab = tree.bab( false, 1e-2);
      auto t2 = chrono::high_resolution_clock::now();
      for_each(x_bab.begin(), x_bab.end(), [](double val) { cout << val << ' '; });
      cout << "\ncx + Q(x) = " << problem.evaluate(x_bab.data()) << '\n';
      cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';


      /*
      DeqForm DEF(env, problem);
      DEF.d_model.set(GRB_IntParam_OutputFlag, 1);
      DEF.solve(300.0);
      //double *x = DEF.d_xVals;
      //for_each(x, x + problem.d_n1, [](double val) { cout << val << ' '; });
      //cout << '\n';
      cout << DEF.d_objVal<< '\n';
      */


      /*
      Benders ben(env, c_env, problem);
      double lpLB = ben.lpSolve();
      double SB = ben.strong_benders();
      cout << "lshaped: " << lpLB << ". SB: " << SB << '\n';

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












