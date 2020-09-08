#include "benders.h"

Benders::Bounds Benders::hybrid_solve(bool lp_cuts, bool sb_cuts, bool zk_cuts, bool strong_cuts, bool affine,
                                      bool force_int, size_t max_rounds, double upper_bound, double tol)
{
  size_t iter = 0;        // number of benders cuts
  size_t round = 0;       // rounds of gmi cuts added

  double LB;
  double UB = GRB_INFINITY;

  bool branch = false;

  while (true)
  {
    Master::Solution sol = d_master.solve();
    if (sol.infeasible)
    {
      cout << "mp infeasible" << endl;
      return Bounds{GRB_INFINITY, GRB_INFINITY};
    }

    vector<double> x = sol.xVals;
    LB = get_lb();

    if (LB > upper_bound)// || LB > UB - 1e-8)
    {
      cout << "LB > upper_bound (LB = " << LB << ")" << endl;
      break;
    }

    bool int_feas = all_of(x.begin(), x.begin() + d_p1, [](double val){ return is_integer(val); });


    if (not int_feas && round < max_rounds)
    {
      cout << "computing gmi cuts\n";
      if (not round_of_cuts(sol, 1e-6))      // at least one cut was added
      {
        ++round;
        cout << "gmi cut added\n";
        continue;
      }
    }

    if (not int_feas)
    {
      cout << "not integer feasible" << endl;
      copy(x.begin(), x.end(), d_xvals);
      branch = true;
      if (force_int)
        break;
    }

    double cx = inner_product(d_problem.d_c.data(), d_problem.d_c.data() + d_n1, x.begin(), 0.0);
    vector<double> vx = d_agg.compute_vwx(x.data());
    double Qx = accumulate(vx.begin(), vx.end(), 0.0) / d_S;
    if (int_feas && cx + Qx < UB)
    {
      copy(x.begin(), x.end(), d_incumbent);
      UB = cx + Qx;
    }
    cout << "LB: " << LB << ". UB: " << UB << endl;

    BendersCut cut;
    if (lp_cuts)
    {
      cout << "computing lp cut\n";
      cut = lpCut(x.data());
      if (not add_cut(cut, sol, tol))
        continue;
    }
    if (sb_cuts)
    {
      //cut = sb_cut(x.data());
      cut = lr_cut(x.data(), vx);
      if (not add_cut(cut, sol, tol))
        continue;
    }
    if (zk_cuts)
    {
      cout << "computing zk cut\n";
      cut = d_pslp.best_zk_cut(sol, d_master, 25, false);
      //cut = d_agg.bac_cut(sol, d_master, tol);
      if (not add_cut(cut, sol, tol))
      {
        cout << "zk cut added\n";
        continue;
      }
    }

    if (strong_cuts)
    {
      cut = d_agg.strong_cut(sol, vx, affine, tol);
      if (not add_cut(cut, sol, tol))
      {
        cout << "strong cut added\n";
        continue;
      }
    }

    copy(x.begin(), x.end(), d_xvals);
    cout << "no improvement possible\n";
    branch = true;
    break;
  }
  //cout << "Number of hybrid cuts: " << iter << '\n';
  cout << "LB: " << LB << ". UB: " << UB << '\n';
  return Bounds { LB, UB, branch };
}





