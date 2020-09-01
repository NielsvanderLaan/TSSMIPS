#include "run.h"

void solve_ri(Data &rand, GRBEnv &env, GRBenv *c_env)
{
  vector<double> eta_star;
  vector<vector<double>> lbda, lbda_plus, gaps, gaps_plus;

  for (size_t instance = 0; instance != 10; ++instance)
  {
    Problem problem(10, 0, 0, 5, 5, 5, 100, rand, env, 0, 0, 0, 5);
    problem.randomInstance(0.1);

    DeqForm DEF(env, problem);
    DEF.solve(300.0);
    double eta = DEF.d_objVal;
    cout << "eta = " << eta;
    eta_star.push_back(eta);

    vector<double> lbda_row, lbda_plus_row, gaps_row, gaps_plus_row;

    Benders ben(env, c_env, problem);
    ben.lpSolve();
    vector<double> alpha(problem.d_m2);
    for (size_t row = 0; row != alpha.size(); ++row)
      alpha[row] = inner_product(problem.d_Tmat[row].begin(), problem.d_Tmat[row].end(), DEF.d_xVals, 0.0);
    ben.lbda(alpha.data(), 1.0);
    cout << ". cxa + Q(xa) = " << problem.evaluate(ben.d_xvals);
    ben.ldSolve(false);
    cout << ". cxa+ + Q(xa+) = " << problem.evaluate(ben.d_xvals) << '\n';



    /*
    for (size_t run = 0; run != 10; ++run)
    {
      vector<double> alpha = rand.unif_real_vec(problem.d_m2, 0, 10);
      Benders copy = ben;
      copy.lbda(alpha.data(), 1.0);
      double cxaQxa = problem.evaluate(copy.d_xvals);
      lbda_row.push_back(cxaQxa);
      gaps_row.push_back((cxaQxa - eta) / eta * 100);

      copy.ldSolve(false);
      double cxaQxa_plus = problem.evaluate(copy.d_xvals);
      lbda_plus_row.push_back(cxaQxa_plus);
      gaps_plus_row.push_back((cxaQxa_plus - eta) / eta * 100);
      cout << (cxaQxa - cxaQxa_plus) / eta * 100 << ' ';
    }
    cout << '\n';
     */
    lbda.push_back(lbda_row);
    lbda_plus.push_back(lbda_plus_row);
    gaps.push_back(gaps_row);
    gaps_plus.push_back(gaps_plus_row);
  }

  /*
  for (auto row : gaps)
  {
    for (auto el : row)
      cout << el << ' ';
    cout << '\n';
  }

  cout << '\n';

  for (auto row : gaps_plus)
  {
    for (auto el : row)
      cout << el << ' ';
    cout << '\n';
  }
  */


}