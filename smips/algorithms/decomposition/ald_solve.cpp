#include "benders.h"

void Benders::ald_solve(double tol, size_t maxRounds)
{
  bool stop = false;
  size_t iter = 0;

  while (not stop)
  {
    cout << "iteration: " << iter << '\n';
    ++iter;
      // solve master problem, and collect x and theta
    Master::Solution sol = d_master.solve();
    // double runtime = d_master.d_model.get(GRB_DoubleAttr_Runtime);
    // cout << "SOLVING MASTER TOOK: " << runtime << " seconds\n"; 
    // bottom line: this becomes expensive (really quickly)

    vector<double> x = sol.xVals;  
    for (size_t var = 0; var != d_n1; ++var)
      cout << x[var] << ' ';
    cout << '\n';
  
    double theta = sol.thetaVal;        
      // derive ald cut
    double beta[d_n1], tau, gamma;
    ald_cut(x.data(), beta, tau, gamma, maxRounds);
    cout << "beta:";
    for (size_t var = 0; var != d_n1; ++var)
      cout << ' ' << beta[var];
    cout << ", tau: " << tau << ", gamma: " << gamma << '\n';

     
      // add ald cut 
    stop = d_master.add_ald_cut(beta, gamma, tau, x.data(), theta, tol);
    
    vector<double> c = d_problem.d_c;
    double cx = 0;
    for (size_t var = 0; var != d_n1; ++var)
      cx += c[var] * x[var];    
    cout << "true Q(x) = " << d_problem.evaluate(x.data()) - cx << "\n\n";
    
    if (stop)
      copy(x.begin(), x.end(), d_xvals);
  }
  
  cout << "Number of ald cuts: " << iter << '\n'; 
}