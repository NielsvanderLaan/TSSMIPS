#include "master.h"

BendersCut Master::fenchel_cut(Solution sol, double tol)
{
  return d_fenchel.fenchel_cut(sol.xVals, sol.thetaVal, tol);
}