#include "benders.h"

double Benders::zk_solve(double tol, size_t maxRounds)
{
  bool stop = false;
  size_t iter = 0;

  while (not stop)
  {
    //cout << "zk iteration: " << iter << '\n';
    ++iter;
      // solve master problem, and collect x and theta
    Master::Solution sol = d_master.solve();

    
    vector<double> x = sol.xVals;  
    
    /*
    cout << "current solution: ";
    for (size_t var = 0; var != d_n1; ++var)
      cout << x[var] << ' ';
    cout << '\n';
    */

    double cx = 0;
    for (size_t var = 0; var != d_n1; ++var)
      cx += d_problem.d_c[var] * x[var];
    double Q = d_problem.evaluate(x.data()) - cx;
    //cout << "Q(x) = " << Q << '\n';

    
    double theta = sol.thetaVal;        
      
      // derive zk cut
    BendersCut cut = d_pslp.best_zk_cut(sol, d_master, maxRounds);

    stop = add_cut(cut, sol, tol);
    
    if (stop)
      copy(x.begin(), x.end(), d_xvals);
    
  }
  
  cout << "Number of zk cuts: " << iter << '\n'; 
  
  double obj;
  GRBgetdblattr(d_master.d_cmodel, "ObjVal", &obj);
  
  return obj;
}





