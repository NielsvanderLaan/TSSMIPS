#include "benders.h"

Benders::Bounds Benders::hybrid_solve(double global_UB, bool affine, double tol)
{
  bool stop = false;
  size_t iter = 0;        // number of benders cuts
  size_t round = 0;       // rounds of gmi cuts added
  size_t count = 0;

  double LB;
  double UB = GRB_INFINITY;
  int max_rounds = 1;
  
  while (not stop)
  {
    Master::Solution sol = d_master.solve();
    if (sol.infeasible)
      return Bounds { GRB_INFINITY, GRB_INFINITY, true };
      
    vector<double> x = sol.xVals;
    LB = get_lb();
    
    if (LB > global_UB - tol)
      break;

    bool int_feas = all_of(x.begin(), x.begin() + d_p1, is_integer);
    if (not int_feas && round < max_rounds)
    {
      if (round_of_cuts(sol, tol))      // no cuts were added
        round = max_rounds;
      else
        ++round;
      continue;
    }

    if (not int_feas)
    {
      copy(x.begin(), x.end(), d_xvals);
      break;
    }

    double cx = inner_product(d_problem.d_c.data(), d_problem.d_c.data() + d_n1, x.begin(), 0.0);
    double Qx = GRB_INFINITY;
    BendersCut cut = d_agg.strong_cut(sol, Qx, affine);

    if (cx + Qx < UB && int_feas)
    {
      copy(x.begin(), x.end(), d_incumbent);
      UB = cx + Qx;
    }
    
    stop = add_cut(cut, sol, tol);
    if (stop)
      copy(x.begin(), x.end(), d_xvals);
    else
      ++iter;
    cout << "LB: " << LB << ". UB: " << UB << '\n';

  }
  cout << "Number of hybrid cuts: " << iter << '\n'; 
  cout << "LB: " << LB << ". UB: " << UB << '\n'; 
 
  return Bounds { LB, UB, false };
}





