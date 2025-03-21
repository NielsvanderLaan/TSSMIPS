#include "zk.h"

bool ZK::solve(double *x, double theta, double rho, Master &master, size_t maxRounds,
               double &nrounds, double &lptime, double &cuttime,
               bool gomory, double tol)
{
  bool stop = false;
  size_t round = 0;

  while (not stop)
  {
        // solve the model by calling optimize(), which also updates d_objval and d_yvals
    auto lp_t1 = chrono::high_resolution_clock::now();
    bool feas = optimize();
    auto lp_t2 = chrono::high_resolution_clock::now();
    lptime += chrono::duration_cast<chrono::milliseconds>(lp_t2 - lp_t1).count()  / 1000.0;

    if (not feas)  // if model is infeasible
      return false;      // return false

    ++round;
    stop = true;        // resets to false if we add a cut (which we attempt to do if solution violates integer constraints)
    if (round >= maxRounds)
      break;
 
    int bhead[d_nConstrs];
    GRBgetBasisHead(d_model, bhead);        // extract basis info   

    size_t nCuts = 0;

    for (size_t row = 0; row != d_nConstrs; ++row)    // loop over rows of simplex tableau 
    {
      int basic_var = bhead[row];            // index of corresponding basic variable
      if (basic_var >= d_p2)                 // check if variable has to be integer
        continue;                            // if not, do not derive a cut

      double yval = d_yvals[basic_var];

      if (is_integer(yval))                  // if variable value is integer,
        continue;                            // then do not derive a cut

      auto t1 = chrono::high_resolution_clock::now();
      Cut cut = gomory ? generate_gmi_cut(master, row, yval, x, theta, rho) : d_cglp.generate_cut(x, rho, d_yvals.data(), basic_var, floor(yval));
      auto t2 = chrono::high_resolution_clock::now();
      cuttime += chrono::duration_cast<chrono::milliseconds>(t2 - t1).count()  / 1000.0;

      if (add_cut(cut, x, rho, tol, d_nConstrs + nCuts, d_nVars + nCuts))    // ret = true iff cut was added (iff cut is proper)
      {
        ++nCuts;
        stop = false;
      }
    }
    d_nConstrs += nCuts;
    d_nVars += nCuts;     // slacks
  }
  nrounds = round;
  return true; // model is feasible
}










