#include "benders.h"

void Benders::lbda(double *alpha, double gomoryTimeLimit, double tol)
{
  d_gomory.setTimeLimit(gomoryTimeLimit);
  
  
  bool stop = false;
  size_t iter = 0;
  
  while (not stop)
  {
    ++iter;
      // solve master problem, and collect x and theta
    Master::Solution sol = d_master.solve(tol);

    vector<double> x = sol.xVals;      
    double theta = sol.thetaVal;        
      // derive cut
    BendersCut cut = lbdaCut(x.data(), alpha);  
    
      // add the cut (conditional on it being violated by the current solution)
    stop = add_cut(cut, sol, tol);  // if no cut was added, then while loop is exited
    
    if (stop)
      copy(x.begin(), x.end(), d_xvals);
  }
  
  //cout << "Number of LBDA iterations: " << iter << '\n';
}











