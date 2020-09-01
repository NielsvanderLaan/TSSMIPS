#include "benders.h"

double Benders::ldSolve(bool affine, double tol)
{
  bool stop = false;
  size_t iter = 0;        // number of benders cuts
  size_t round = 0;       // rounds of gmi cuts added
  size_t count = 0;

  double LB;

  while (not stop)
  {
    //cout << "iter = " << iter << '\n';
    Master::Solution sol = d_master.solve();
    vector<double> x = sol.xVals;
    LB = get_lb();
    //cout << "LB = " << LB << '\n';

    vector<double> vx = d_agg.compute_vwx(x.data());
    double Qx = accumulate(vx.begin(), vx.end(), 0.0) / d_S;
    BendersCut cut;
    if (affine)
      cut = lr_cut(x.data(), vx);
    else
      cut = d_agg.strong_cut(sol, vx, affine, tol);
    //cout << "Q(x) = " << Qx << '\n';
    //cout << "unscaled cut value = " << -inner_product(cut.d_beta.begin(), cut.d_beta.end(), x.data(), -cut.d_alpha) - cut.d_tau * sol.thetaVal << '\n';
    //cout << "tau = " << cut.d_tau << '\n';
    //cout << "scaled cut value = " << -inner_product(cut.d_beta.begin(), cut.d_beta.end(), x.data(), -cut.d_alpha) / (1 + cut.d_tau) << '\n';

    stop = add_cut(cut, sol, tol);
    if (stop)     // adding strong cut
      copy(x.begin(), x.end(), d_xvals);
    else
      ++iter;


  }

  //cout << "number of cuts: " << iter << '\n';


  return LB;
}