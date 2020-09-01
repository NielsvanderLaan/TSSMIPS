#include "zk.h"

bool ZK::solve(double *x, double theta, Master &master, size_t maxRounds, bool gomory, bool zk, double tol)
{
  bool stop = false;
  size_t round = 0;
      
  while (not stop)
  {
        // solve the model by calling optimize(), which also updates d_objval and d_yvals
    if (not optimize())  // if model is infeasible
      return false;      // return false

    ++round;
    stop = true;        // resets to false if we add a cut (which we attempt to do if solution violates integer constraints)
    if (round >= maxRounds)
      break;
 
    int bhead[d_nConstrs];
    GRBgetBasisHead(d_model, bhead);        // extract basis info   

    
    /*
    for (size_t var = 0; var != d_n2; ++var)
       cout << d_yvals[var] << ' ';
    cout << endl;
    double obj;
    GRBgetdblattr(d_model, "ObjVal", &obj);
    cout << "objective = " << obj << '\n';
    */


    size_t nCuts = 0;
    for (size_t row = 0; row != d_nConstrs; ++row)    // loop over rows of simplex tableau 
    {
      int basic_var = bhead[row];            // index of corresponding basic variable
      if (basic_var >= d_p2)                 // check if variable has to be integer
        continue;                            // if not, do not derive a cut

      double yval = d_yvals[basic_var];

      if (is_integer(yval))                  // if variable value is integer,
        continue;                            // then do not derive a cut

      Cut cut = gomory ? generate_gmi_cut(master, row, yval, x, zk) : d_cglp.generate_cut(x, theta, d_yvals.data(), basic_var, floor(yval));


      if (add_cut(cut, x, theta, tol, d_nConstrs + nCuts))    // ret = true iff cut was added (iff cut is proper)
      {
        ++nCuts;
        stop = false;
      }

    }
    d_nConstrs += nCuts;
    d_nVars += nCuts;     // slacks
  }
  return true; // model is feasible
}










