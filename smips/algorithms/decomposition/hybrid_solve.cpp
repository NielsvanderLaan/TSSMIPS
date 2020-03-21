#include "benders.h"

Benders::Bounds Benders::hybrid_solve(double global_UB, bool affine, double tol, int max_iter)
{
  bool stop = false;
  size_t iter = 0;        // number of benders cuts
  size_t count = 0;       // number of cuts since integer feasible solution was encountered 
  size_t round = 0;       // rounds of gmi cuts added

  double LB;
  double UB = 1e20;
  int max_rounds = 3;
  
  while (not stop)
  {
    Master::Solution sol = d_master.solve();
    if (sol.infeasible)
      return Bounds { 1e20, 1e20, true };
      
    vector<double> x = sol.xVals;  
    LB = get_lb();
    
    if (LB > global_UB - tol)
      break;
  
    bool int_feas = all_of(x.begin(), x.begin() + d_p1, is_integer);
    if (not int_feas && round < max_rounds)
    {
      if (round_of_cuts(sol, tol))
        round = max_rounds;
      else
        ++round;
      continue;
    }
    
    count = int_feas ? 0 : count + 1;    // reset counter if current solution is integer
    if (count > max_iter && not int_feas)
    {
      copy(x.begin(), x.end(), d_xvals);
      break;
    }
    
    double cx = inner_product(d_problem.d_c.data(), d_problem.d_c.data() + d_n1, x.begin(), 0.0);
    double Qx = 0.0;
    BendersCut cut = d_agg.strong_cut(sol, Qx, affine);

    if (cx + Qx < UB && int_feas)      
    {
      copy(x.begin(), x.end(), d_incumbent);
      UB = cx + Qx;
    }
    
    stop = add_cut(cut, sol, tol);
    if (stop)
      copy(x.begin(), x.end(), d_xvals);

    ++iter;
    cout << "LB: " << LB << ". UB: " << UB << '\n';

  }
  cout << "Number of hybrid cuts: " << iter << '\n'; 
  cout << "LB: " << LB << ". UB: " << UB << '\n'; 
 
  return Bounds { LB, UB, false };
}





