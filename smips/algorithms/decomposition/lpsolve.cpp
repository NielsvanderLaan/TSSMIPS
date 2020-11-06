#include "benders.h"

double Benders::lpSolve(double tol)
{
  bool stop = false;
  size_t iter = 0;
  
  while (not stop)
  {
    Master::Solution sol = d_master.solve(tol);
  
    BendersCut cut = d_agg.lp_cut( sol.xVals);
    stop = add_cut(cut, sol, tol);

    if (stop)
      copy( sol.xVals.begin(),  sol.xVals.end(), d_xvals);
    else
      ++iter;
  }
  
  cout << "lp cuts: " << iter << '\n';
  double obj;
  GRBgetdblattr(d_master.d_cmodel, "ObjVal", &obj);
  return obj;
}











