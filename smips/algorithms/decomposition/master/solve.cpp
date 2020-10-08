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
  if (status != 2)
  {
    cerr << "master problem status: " << status << '\n';
    exit(status);
  }

  if (d_zk_safe)
  {
    double violation, resid;
    GRBgetdblattr(d_cmodel, "ConstrVio", &violation);
    GRBgetdblattr(d_cmodel, "ConstrResidual", &resid);

    if (violation + resid > 1e-4)
    {
      cout << "master violation (before) = " << violation << ", resid = " << resid << '\n';

      chg_mp_tol(true);
      GRBoptimize(d_cmodel);
      chg_mp_tol(false);

      GRBgetdblattr(d_cmodel, "ConstrVio", &violation);
      GRBgetdblattr(d_cmodel, "ConstrResidual", &resid);
      cout << "master violation (after) = " << violation << ", resid = " << resid << '\n';

      GRBgetintattr(d_cmodel, "Status", &status);
      if (status == 3 || status == 4)
        return Solution{ vector<double>(0), -1, true };
      if (status == 5)
      {
        cout << "master problem status = 5\n";
        return Solution{ vector<double>(0), -1, true };
      }
      if (status != 2)
      {
        cerr << "master problem status: " << status << '\n';
        exit(status);
      }
    }
  }

  vector<double> x(d_n1);
  GRBgetdblattrarray(d_cmodel, "X", 1, d_n1, x.data());
  double theta;
  GRBgetdblattrelement(d_cmodel, "X", 0, &theta);

  return Solution{ x, theta + d_L, false };
}