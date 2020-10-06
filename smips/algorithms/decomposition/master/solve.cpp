#include "master.h"

#include <string>

Master::Solution Master::solve(double tol)
{
  GRBoptimize(d_cmodel);

  int status;
  GRBgetintattr(d_cmodel, "Status", &status);

  if (status == 3 || status == 4)      // model is infeasible
    return Solution{ vector<double>(0), -1, true };


  if (d_zk_safe)
  {
    double violation, resid;
    GRBgetdblattr(d_cmodel, "ConstrVio", &violation);
    GRBgetdblattr(d_cmodel, "ConstrResidual", &resid);

    if (violation + resid > 1e-4)
    {
      cout << "master violation (before) = " << violation << ", resid = " << resid << '\n';

      GRBreset(d_cmodel, 0);
      GRBenv *env = GRBgetenv(d_cmodel);
      GRBsetintparam(env, "NumericFocus", 3);
      GRBsetintparam(env, "ScaleFlag", 0);
      GRBsetintparam(env, "Method", 1);
      GRBsetintparam(env, "PreSolve", 1);
      GRBoptimize(d_cmodel);
      GRBsetintparam(env, "NumericFocus", 0);
      GRBsetintparam(env, "ScaleFlag", -1);
      GRBsetintparam(env, "Method", -1);
      GRBsetintparam(env, "PreSolve", -1);

      GRBgetdblattr(d_cmodel, "ConstrVio", &violation);
      GRBgetdblattr(d_cmodel, "ConstrResidual", &resid);
      cout << "master violation (after) = " << violation << ", resid = " << resid << '\n';
    }
  }
  //GRBwrite(d_cmodel, "master.lp");



  vector<double> x(d_n1);
  GRBgetdblattrarray(d_cmodel, "X", 1, d_n1, x.data());   
  
  double theta;
  GRBgetdblattrelement(d_cmodel, "X", 0, &theta);
  return Solution{ x, theta + d_L, false };
}