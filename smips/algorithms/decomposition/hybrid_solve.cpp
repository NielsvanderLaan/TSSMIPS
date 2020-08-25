#include "benders.h"

Benders::Bounds Benders::hybrid_solve(double upper_bound, bool affine, bool lp_cuts, double tol)
{
  size_t iter = 0;        // number of benders cuts
  size_t round = 0;       // rounds of gmi cuts added

  double LB;
  double UB = GRB_INFINITY;
  int max_rounds = 25;

  bool branch = false;

  while (true)
  {
    Master::Solution sol = d_master.solve();
    if (sol.infeasible)
    {
      //cout << "mp infeasible\n";
      return Bounds{GRB_INFINITY, GRB_INFINITY};
    }

    vector<double> x = sol.xVals;
    LB = get_lb();

    if (LB > upper_bound)
    {
      //cout << "LB > upper_bound (LB = " << LB << ")\n";
      break;
    }

    bool int_feas = all_of(x.begin(), x.begin() + d_p1, [](double val) { return is_integer(val); });

    if (not int_feas && round < max_rounds)
    {
      if (not round_of_cuts(sol, tol))      // at least one cut was added
      {
        ++round;
        continue;
      }
    }

    if (not int_feas)
    {
      copy(x.begin(), x.end(), d_xvals);
      //cout << "not integer feasible\n";
      branch = true;
      break;
    }

    double cx = inner_product(d_problem.d_c.data(), d_problem.d_c.data() + d_n1, x.begin(), 0.0);
    vector<double> vx = d_agg.compute_vwx(x.data());
    double Qx = accumulate(vx.begin(), vx.end(), 0.0) / d_S;
    if (cx + Qx < UB)     // this point is only reached if x is integer feasible
    {
      copy(x.begin(), x.end(), d_incumbent);
      UB = cx + Qx;
    }
    //cout << "LB: " << LB << ". UB: " << UB << '\n';

    BendersCut cut;
    /*
    if (lp_cuts)
    {
      cut = lpCut(x.data());
      if (not add_cut(cut, sol, tol))
        continue;
    }
     */

    //cut = d_agg.bac_cut(sol, d_master, tol, 10);

    //cut = d_pslp.best_zk_cut(sol, d_master, 10, false);
    cut = sb_cut(x.data());
    if (not add_cut(cut, sol, tol))
      continue;
    else
    {
      copy(x.begin(), x.end(), d_xvals);
      //cout << "no improvement possible\n";
      branch = true;
      break;
    }


    cut = d_agg.strong_cut(sol, vx, affine, tol);

    if (add_cut(cut, sol, tol))     // adding strong cut
    {
      copy(x.begin(), x.end(), d_xvals);
      cout << "no improvement possible\n";
      branch = true;
      break;
    } else
      ++iter;
  }
  //cout << "Number of hybrid cuts: " << iter << '\n';
  //cout << "LB: " << LB << ". UB: " << UB << '\n';

  return Bounds{LB, UB, branch};
}





