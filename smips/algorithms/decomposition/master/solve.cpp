#include "master.h"

#include <string>

Master::Solution Master::solve()
{
  GRBoptimize(d_cmodel);       
  
  int status;
  GRBgetintattr(d_cmodel, "Status", &status);

  if (status == 3 || status == 4)      // model is infeasible
    return Solution{ vector<double>(0), -1, true };

  vector<double> x(d_n1);
  GRBgetdblattrarray(d_cmodel, "X", 1, d_n1, x.data());   
  
  double theta;
  GRBgetdblattrelement(d_cmodel, "X", 0, &theta);
  return Solution{ x, theta, false };
}