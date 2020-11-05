#include "zk.h"

BendersCut ZK::subgradient()
{
  double lambda[d_nConstrs];
  GRBgetdblattrarray(d_model, "Pi", 0, d_nConstrs, lambda);
  
  double alpha = 0; 
  vector<double> beta(d_n1, 0.0);
  double tau = 0;
  
  for (size_t con = 0; con != d_nConstrs; ++con)
  { 
    alpha += lambda[con] * d_omega[con];
    tau += lambda[con] * d_tau[con]; 
    
    for (size_t var = 0; var != d_n1; ++var)
      beta[var] += lambda[con] * d_Tmat[con][var];  
  }

  return BendersCut{ alpha, beta, tau };
}