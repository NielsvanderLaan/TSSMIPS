#include "benders.h"

double Benders::strong_benders(double tol)
{
  bool stop = false;
  
  size_t iter = 0;
  
  while (not stop)
  {
    ++iter;
      // solve master problem, and collect x and theta
    Master::Solution sol = d_master.solve();
    vector<double> x = sol.xVals;      
      // derive cut
    BendersCut cut = sb_cut(x.data());  
      // add the cut (conditional on it being violated by the current solution)
    stop = add_cut(cut, sol, tol);  // if no cut was added, then while loop is exited
    
    if (stop)
      copy(x.begin(), x.end(), d_xvals);  
  }
  
  cout << "Number of strong benders cuts: " << iter << '\n';
  double obj;
  GRBgetdblattr(d_master.d_cmodel, "ObjVal", &obj);
  return obj;
}











