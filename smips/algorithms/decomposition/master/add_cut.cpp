#include "master.h"


bool Master::add_cut(BendersCut cut, Solution sol, double tol)
{
      // adding the constraint: (1 + tau) theta + beta^T x >= alpha  
  double alpha_betax = -inner_product(cut.d_beta.begin(), cut.d_beta.end(), sol.xVals.begin(), -cut.d_alpha);    // alpha - beta^T x

  if (alpha_betax > (1 + cut.d_tau) * sol.thetaVal + tol) // then add cut and return false
  {
    add_cut(cut);
    return false;
  }
  return true; // betaxgamma >= theta, no cut added 
}