#include "master.h"

#include <string>

Master::Solution Master::solve(double tol)
{
  GRBoptimize(d_cmodel);

  int status;
  GRBgetintattr(d_cmodel, "Status", &status);

  if (status == 3 || status == 4)
    return Solution{ vector<double>(0), -1, true };

  if (status == 5)
  {
    cout << "master problem status = 5\n";
    return Solution{ vector<double>(0), -1, true };
  }

  double vio = violation();
  if (vio > 1e-4 || status != 2)
  {
    cout << "master violation (before) = " << vio << ". status: " << status << '\n';

    chg_mp_tol(true);
    GRBoptimize(d_cmodel);
    chg_mp_tol(false);

    vio = violation();
    GRBgetintattr(d_cmodel, "Status", &status);
    cout << "master violation (after) = " << vio << ". status: " << status << '\n';

    if (status == 3)
    {
      GRBwrite(d_cmodel, "master.lp");
      chg_mp_tol(true);
      GRBsetintparam(GRBgetenv(d_cmodel), "OutputFlag", 1);
      GRBoptimize(d_cmodel);
      exit(3);
    }


    if (status == 3 || status == 4 || status == 5)
      return Solution{ vector<double>(0), -1, true };
    if (status != 2)
    {
      cout << "master problem status: " << status << '\n';
      exit(status);
    }
  }

  vector<double> x(d_n1);
  GRBgetdblattrarray(d_cmodel, "X", 1, d_n1, x.data());
  double theta;
  GRBgetdblattrelement(d_cmodel, "X", 0, &theta);

  return Solution{ x, theta + d_L, false };
}