#include "benders.h"

void Benders::lbda(vector<double> &alpha, double gomoryTimeLimit, double tol)
{
  d_gomory.setTimeLimit(gomoryTimeLimit);

  bool stop = false;
  size_t iter = 0;
  
  while (not stop)
  {
    ++iter;
    Master::Solution sol = d_master.solve(tol);

    BendersCut cut = lbdaCut(sol.xVals, alpha);

    stop = add_cut(cut, sol, tol);
    
    if (stop)
      copy(sol.xVals.begin(), sol.xVals.end(), d_xvals);
  }
  
  //cout << "Number of LBDA iterations: " << iter << '\n';
}











