#include "cgmip.h"

BendersCut CGMip::get_candidate()
{  
  double *beta_ptr = d_mp.get(GRB_DoubleAttr_X, d_beta.data(), d_beta.size());
  vector<double> beta(beta_ptr, beta_ptr + d_beta.size());
  delete[] beta_ptr;     
  
                             // we set tau = max(d_tau, 0) to prevent numerical issues                  
  return BendersCut{ d_alpha.get(GRB_DoubleAttr_X), beta, max(d_tau.get(GRB_DoubleAttr_X), 0.0) }; 
}