#include "run.h"

void lbda_scheme(Problem &problem, GRBEnv &env, GRBenv *c_env, size_t nIter)
{
  Benders lshaped(env, c_env, problem, false);
  lshaped.lpSolve();
  Benders ben = lshaped;
  vector<double> alpha(problem.d_m2);
  ben.lbda(alpha.data(), 1.0);
  double lbda = problem.evaluate(ben.d_xvals);
  ben.ldSolve(false, 0.01);
  double lbda_plus = problem.evaluate(ben.d_xvals);

  double best_lbda = lbda;
  double best_lbda_plus = lbda_plus;

  cout << "cx_a + Q(x_a) = " << lbda << ".\ncx_a^+ + Q(x_a^+) = " << lbda_plus << '\n';


  cout << "best lbda = " << best_lbda << ".\nBest strong lbda = " << best_lbda_plus << '\n';


  for (size_t iter = 0; iter != nIter; ++iter)
  {
    for (size_t idx = 0; idx != alpha.size(); ++idx)
      alpha[idx] = inner_product(problem.d_Tmat[idx].begin(), problem.d_Tmat[idx].end(), ben.d_xvals, 0.0);
    cout << "alpha = ";
    for_each(alpha.begin(), alpha.end(), [](double val){cout << val << ' ';});
    cout << endl;


    Benders ben = lshaped;

    ben.lbda(alpha.data(), 1.0);
    double lbda = problem.evaluate(ben.d_xvals);
    ben.ldSolve(false, 0.01);
    double lbda_plus = problem.evaluate(ben.d_xvals);

    best_lbda = min(best_lbda, lbda);
    best_lbda_plus = min(best_lbda_plus, lbda_plus);

    cout << "cx_a + Q(x_a) = " << lbda << ".\ncx_a^+ + Q(x_a^+) = " << lbda_plus << '\n';
  }


}