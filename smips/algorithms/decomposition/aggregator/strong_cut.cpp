#include "aggregator.h"

BendersCut Aggregator::strong_cut(Master::Solution sol, double &Qx, double affine, double tol)
{ 
  double rho = sol.thetaVal;
  double *x = sol.xVals.data();
  double cRho = 1;

  BendersCut cut;
  
  bool first_time = true;
  while (cRho > tol)
  {
    cRho = -rho;
    cut = BendersCut{ 0, vector<double>(d_n1), 0};

    for (size_t s = 0; s != d_cgmips.size(); ++s)
    {
      double prob = d_probs[s];
      
      double vwx = -1;
      if (first_time)
      {
        vwx = compute_vwx(x, s);
        Qx += prob * vwx;  
      }
      cut += d_cgmips[s].generate_cut(x, rho, first_time, vwx, affine) * prob; 
      cRho -= prob * d_cgmips[s].mp_val();
    }  
    rho += cRho / (1 + cut.d_tau);    
    first_time = false;
  } 

  return cut;
}


